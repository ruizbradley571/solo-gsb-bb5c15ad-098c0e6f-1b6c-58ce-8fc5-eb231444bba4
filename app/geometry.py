"""二维几何计算与斯涅尔定律折射。

约定：
- 所有方向向量、法线在内部均归一化。
- 界面为一条直线：由线上一点 point 与指向下一介质的法线 normal 定义。
- 光线沿当前方向前进，只接受前向（t > 0）交点。
"""
from __future__ import annotations

import math
from dataclasses import dataclass

EPS = 1e-12
T_MIN = 1e-9  # 交点参数 t 必须大于该值，避免自相交

AIM_DEV_TOL = 1e-6  # 目标偏差容差：|偏差| 不大于该值即视为命中
AIM_ANGLE_TOL_DEG = 1e-7  # 角度区间宽容差（度）：二分停止条件，兼作解去重阈值
AIM_MAX_ITER = 60  # 单个异号区间的最大二分迭代次数

FOCUS_MIN_HITS = 2  # 扫描/细化位置有效所需的最少命中光线数
FOCUS_MAX_ROUNDS = 24  # 单个低谷区间黄金分割细化的最多轮数
# 黄金分割内点位置：左内点系数 1/φ²≈0.381966，右内点系数 1/φ≈0.618034
_GOLDEN_TAU = (math.sqrt(5.0) - 1.0) / 2.0  # 1/φ ≈ 0.618034
_GOLDEN_1MTAU = 1.0 - _GOLDEN_TAU  # 1/φ² ≈ 0.381966


def _dot(a: tuple[float, float], b: tuple[float, float]) -> float:
    return a[0] * b[0] + a[1] * b[1]


def normalize(v: tuple[float, float]) -> tuple[float, float]:
    n = math.hypot(v[0], v[1])
    if n < EPS:
        raise ValueError("zero-length vector")
    return (v[0] / n, v[1] / n)


@dataclass
class SegmentResult:
    """一段光线的计算结果。hit_point 为 None 表示光线出射（无后续界面）。"""

    origin: tuple[float, float]
    direction: tuple[float, float]
    n_medium: float
    hit_point: tuple[float, float] | None
    event: str | None  # refraction / total_internal_reflection / no_intersection / None
    n_after: float | None


def _intersect_line(
    origin: tuple[float, float],
    direction: tuple[float, float],
    point: tuple[float, float],
    normal: tuple[float, float],
) -> tuple[float, float] | None:
    """光线与界面直线求交，返回前向交点；平行或交点在身后则返回 None。"""
    denom = _dot(direction, normal)
    if abs(denom) < EPS:
        return None  # 光线与界面平行
    t = _dot((point[0] - origin[0], point[1] - origin[1]), normal) / denom
    if t <= T_MIN:
        return None  # 交点不在前方
    return (origin[0] + t * direction[0], origin[1] + t * direction[1])


def _refract(
    direction: tuple[float, float],
    normal: tuple[float, float],
    n1: float,
    n2: float,
) -> tuple[float, float] | None:
    """按斯涅尔定律计算折射方向；发生全反射时返回 None。

    若光线从法线背面入射，则翻转法线再计算（折射率仍按给定 n1 -> n2）。
    """
    cos_i = _dot(direction, normal)
    if cos_i < 0.0:
        normal = (-normal[0], -normal[1])
        cos_i = -cos_i
    eta = n1 / n2
    sin2_t = eta * eta * max(0.0, 1.0 - cos_i * cos_i)
    if sin2_t > 1.0:
        return None  # 全反射
    cos_t = math.sqrt(1.0 - sin2_t)
    k = cos_t - eta * cos_i
    return (eta * direction[0] + k * normal[0], eta * direction[1] + k * normal[1])


def fresnel_transmittance(
    direction: tuple[float, float],
    normal: tuple[float, float],
    n1: float,
    n2: float,
) -> tuple[float, float, float, float] | None:
    """单界面菲涅耳 s/p 功率透过率（忽略吸收，T = 1 - R）。

    以光线方向与界面法线求入射角余弦 cos_i（光线从法线背面入射时翻转法线），
    再由斯涅尔定律 sin_t = (n1/n2) sin_i 求折射角余弦 cos_t：
        Rs = [(n1 cos_i - n2 cos_t)/(n1 cos_i + n2 cos_t)]²
        Rp = [(n2 cos_i - n1 cos_t)/(n2 cos_i + n1 cos_t)]²
    发生全反射（sin²t > 1）时返回 None。
    返回 (cos_i, cos_t, 1 - Rs, 1 - Rp)。
    """
    nrm = normalize(normal)
    cos_i = _dot(direction, nrm)
    if cos_i < 0.0:
        nrm = (-nrm[0], -nrm[1])
        cos_i = -cos_i
    eta = n1 / n2
    sin2_t = eta * eta * max(0.0, 1.0 - cos_i * cos_i)
    if sin2_t > 1.0:
        return None  # 全反射
    cos_t = math.sqrt(max(0.0, 1.0 - sin2_t))
    r_s = (n1 * cos_i - n2 * cos_t) / (n1 * cos_i + n2 * cos_t)
    r_p = (n2 * cos_i - n1 * cos_t) / (n2 * cos_i + n1 * cos_t)
    return cos_i, cos_t, 1.0 - r_s * r_s, 1.0 - r_p * r_p


def ray_transmittance(
    segments: list[SegmentResult],
    status: str,
    interfaces: list[tuple[tuple[float, float], tuple[float, float], float]],
) -> tuple[float, float]:
    """对一条已追迹光线逐界面累乘 s/p 功率透过率（各界面 T = 1 - R）。

    复用 trace_ray 的结果：仅对 status=completed 的光线，按界面对应顺序取
    各折射段的入射方向与界面法线计算菲涅耳反射率；全反射或中途漏过界面
    （status 非 completed）时整条光线透过率记 0。无界面时透过率为 (1, 1)。
    """
    if status != "completed":
        return 0.0, 0.0
    t_s, t_p = 1.0, 1.0
    for seg, (_point, normal, n_next) in zip(segments, interfaces):
        fresnel = fresnel_transmittance(seg.direction, normal, seg.n_medium, n_next)
        if fresnel is None:
            return 0.0, 0.0  # 全反射：零功率
        t_s *= fresnel[2]
        t_p *= fresnel[3]
    return t_s, t_p


def sample_directions(
    center: tuple[float, float],
    half_angle_deg: float,
    sample_count: int,
) -> list[tuple[float, tuple[float, float]]]:
    """以中心方向为 0°，在 [-a, +a]（度，含端点）等间距生成 sample_count 个方向。

    返回 [(emission_angle_deg, normalized_direction), ...]，角度正向为数学上的
    逆时针方向，端点精确取 -a 与 +a。
    """
    base = normalize(center)
    theta0 = math.atan2(base[1], base[0])
    a = float(half_angle_deg)
    last = sample_count - 1
    result: list[tuple[float, tuple[float, float]]] = []
    for i in range(sample_count):
        if i == 0:
            angle = -a
        elif i == last:
            angle = a
        else:
            angle = -a + 2.0 * a * i / last
        theta = theta0 + math.radians(angle)
        result.append((angle, (math.cos(theta), math.sin(theta))))
    return result


def observation_intersection(
    segments: list[SegmentResult],
    point: tuple[float, float],
    normal: tuple[float, float],
) -> tuple[float, float] | None:
    """取完成追迹后的出射段（最后一段），求其与观察线的前向交点。

    平行（denom≈0）或交点在光线反方向（t<=0）时返回 None。
    """
    final = segments[-1]
    return _intersect_line(final.origin, final.direction, point, normalize(normal))


def centroid_and_rms(
    points: list[tuple[float, float]],
) -> tuple[tuple[float, float], float] | None:
    """点集质心与 RMS 半径（各点到质心距离平方均值的开方）；空集返回 None。"""
    count = len(points)
    if count == 0:
        return None
    cx = sum(p[0] for p in points) / count
    cy = sum(p[1] for p in points) / count
    mean_sq = sum((p[0] - cx) ** 2 + (p[1] - cy) ** 2 for p in points) / count
    return (cx, cy), math.sqrt(mean_sq)


def weighted_centroid(
    points: list[tuple[tuple[float, float], float]],
) -> tuple[float, float] | None:
    """功率加权质心：sum(w·p)/sum(w)；权重总和为 0（或空集）时返回 None。"""
    total = sum(w for _p, w in points)
    if total <= EPS:
        return None
    cx = sum(w * p[0] for p, w in points) / total
    cy = sum(w * p[1] for p, w in points) / total
    return (cx, cy)


def trace_ray(
    start: tuple[float, float],
    direction: tuple[float, float],
    n_start: float,
    interfaces: list[tuple[tuple[float, float], tuple[float, float], float]],
) -> tuple[list[SegmentResult], str, str | None]:
    """依次对每个界面做前向求交并折射。

    interfaces: 按顺序排列的 (point, normal, n_next)。
    返回 (segments, status, reason)，status ∈
    completed / total_internal_reflection / no_intersection。
    """
    origin = (float(start[0]), float(start[1]))
    d = normalize(direction)
    n_current = float(n_start)
    segments: list[SegmentResult] = []

    for idx, (point, normal, n_next) in enumerate(interfaces):
        nrm = normalize(normal)
        hit = _intersect_line(origin, d, point, nrm)
        if hit is None:
            reason = (
                f"光线与第 {idx} 个界面平行，或交点位于光线反方向，追迹停止"
            )
            segments.append(SegmentResult(origin, d, n_current, None, "no_intersection", None))
            return segments, "no_intersection", reason

        new_dir = _refract(d, nrm, n_current, n_next)
        if new_dir is None:
            reason = (
                f"在第 {idx} 个界面发生全反射：n1={n_current} -> n2={n_next}，"
                f"入射角超过临界角，追迹停止"
            )
            segments.append(
                SegmentResult(origin, d, n_current, hit, "total_internal_reflection", n_next)
            )
            return segments, "total_internal_reflection", reason

        segments.append(SegmentResult(origin, d, n_current, hit, "refraction", n_next))
        origin = hit
        d = new_dir
        n_current = n_next

    # 穿过所有界面，最后一段出射
    segments.append(SegmentResult(origin, d, n_current, None, None, None))
    return segments, "completed", None


def _rotate_normal(
    normal: tuple[float, float], angle_deg: float
) -> tuple[float, float]:
    """将单位法线绕原点旋转 angle_deg 度（逆时针为正）。"""
    theta = math.radians(angle_deg)
    c, s = math.cos(theta), math.sin(theta)
    return (c * normal[0] - s * normal[1], s * normal[0] + c * normal[1])


def alignment_sensitivities(
    start: tuple[float, float],
    direction: tuple[float, float],
    n_start: float,
    interfaces: list[tuple[tuple[float, float], tuple[float, float], float]],
    obs_point: tuple[float, float],
    obs_normal: tuple[float, float],
    perturbation_deg: float,
) -> dict:
    """基准追迹 + 逐界面法线 ±扰动，按中心差分求观察线偏移对法线角度的导数。

    观察线法线归一化后逆时针旋转 90° 得单位切向量 t；偏移量 = 出射段与观察线
    前向交点相对 obs_point 沿 t 的有符号距离。

    对每个界面只将其法线旋转 ±perturbation_deg（界面点与其他参数不变）后重新
    追迹，导数 = (s(+δ) - s(-δ)) / (2δ)，单位为长度/度。基准或任一扰动发生
    全反射、漏过界面或无法前向命中观察线时，该次偏移为 None、对应导数为 None，
    不跨失效分支求差分。该函数为即时计算，不写入数据库。

    返回 dict：
      tangent: 观察线单位切向量；
      baseline: (status, reason, point|None, offset|None)；
      interfaces: 每个界面 [plus, minus, derivative|None]，其中
                  plus/minus = (status, reason, point|None, offset|None, perturbed_normal)。
    """
    nrm = normalize(obs_normal)
    tangent = (-nrm[1], nrm[0])  # 法线逆时针旋转 90° 得单位切向量
    normals = [normalize(normal) for _point, normal, _n in interfaces]

    def run(
        itf: list[tuple[tuple[float, float], tuple[float, float], float]],
    ) -> tuple[str, str | None, tuple[float, float] | None, float | None]:
        """追迹一次并求观察线偏移；全反射/漏过界面/无前向交点时偏移为 None。"""
        segments, status, reason = trace_ray(start, direction, n_start, itf)
        if status != "completed":
            return status, reason, None, None
        point = observation_intersection(segments, obs_point, obs_normal)
        if point is None:
            return (
                "no_observation_intersection",
                "出射段与观察线平行，或交点位于光线反方向",
                None,
                None,
            )
        offset = _dot((point[0] - obs_point[0], point[1] - obs_point[1]), tangent)
        return "completed", None, point, offset

    baseline = run(interfaces)

    results: list[
        tuple[
            tuple[str, str | None, tuple[float, float] | None, float | None, tuple],
            tuple[str, str | None, tuple[float, float] | None, float | None, tuple],
            float | None,
        ]
    ] = []
    for idx, ((point, _normal, n_next), base_nrm) in enumerate(
        zip(interfaces, normals)
    ):
        plus_nrm = _rotate_normal(base_nrm, perturbation_deg)
        minus_nrm = _rotate_normal(base_nrm, -perturbation_deg)
        itf_plus = [
            (pt, plus_nrm if i == idx else nm, nn)
            for i, (pt, nm, nn) in enumerate(interfaces)
        ]
        itf_minus = [
            (pt, minus_nrm if i == idx else nm, nn)
            for i, (pt, nm, nn) in enumerate(interfaces)
        ]
        r_plus = run(itf_plus) + (plus_nrm,)
        r_minus = run(itf_minus) + (minus_nrm,)
        derivative: float | None = None
        if r_plus[3] is not None and r_minus[3] is not None:
            derivative = (r_plus[3] - r_minus[3]) / (2.0 * perturbation_deg)
        results.append((r_plus, r_minus, derivative))

    return {"tangent": tangent, "baseline": baseline, "interfaces": results}


def aim_solutions(
    start: tuple[float, float],
    center_direction: tuple[float, float],
    half_angle_deg: float,
    sample_count: int,
    n_start: float,
    interfaces: list[tuple[tuple[float, float], tuple[float, float], float]],
    obs_point: tuple[float, float],
    obs_normal: tuple[float, float],
    target_offset: float,
) -> list[tuple[float, tuple[float, float], tuple[float, float], float]]:
    """求解使出射光线命中观察线上目标位置的发射角（扫描 + 异号区间二分）。

    目标位置 = obs_point + target_offset · t，其中 t 为观察线单位切向量
    （法线逆时针旋转 90°）；目标偏差 = 观察线交点相对目标位置沿 t 的有符号距离。

    流程：
    - 在 [-half_angle_deg, +half_angle_deg]（度，含端点）内等间距取 sample_count
      个发射角逐一求偏差，|偏差| ≤ 1e-6 的采样点直接记为解；
    - 仅对相邻有效采样（追迹完成且与观察线有前向交点）偏差异号的区间做二分，
      二分过程中遇无效追迹或无前向交点即丢弃该区间，不跨无效采样连根；
    - 二分停止条件：|偏差| ≤ 1e-6 或角度区间宽 ≤ 1e-7 度，最多迭代 60 次。

    返回 [(angle_deg, exit_direction, observation_point, deviation), ...]，
    按角度升序、按 1e-7 度去重（保留 |偏差| 较小者）；无解时为空列表。
    """
    nrm = normalize(obs_normal)
    tangent = (-nrm[1], nrm[0])  # 法线逆时针旋转 90° 得单位切向量
    base = normalize(center_direction)
    theta0 = math.atan2(base[1], base[0])

    def evaluate(
        angle_deg: float,
    ) -> tuple[float, tuple[float, float], tuple[float, float]] | None:
        """指定发射角的目标偏差、出射方向与观察线交点；无效时返回 None。"""
        theta = theta0 + math.radians(angle_deg)
        direction = (math.cos(theta), math.sin(theta))
        segments, status, _ = trace_ray(start, direction, n_start, interfaces)
        if status != "completed":
            return None
        point = observation_intersection(segments, obs_point, obs_normal)
        if point is None:
            return None
        deviation = (
            _dot((point[0] - obs_point[0], point[1] - obs_point[1]), tangent)
            - target_offset
        )
        return deviation, segments[-1].direction, point

    solutions: list[
        tuple[float, tuple[float, float], tuple[float, float], float]
    ] = []

    # 扫描：等间距采样，记录每个采样点的求值结果，直接命中者先记为解
    sampled: list[
        tuple[float, tuple[float, tuple[float, float], tuple[float, float]] | None]
    ] = []
    for angle, _ in sample_directions(center_direction, half_angle_deg, sample_count):
        result = evaluate(angle)
        sampled.append((angle, result))
        if result is not None and abs(result[0]) <= AIM_DEV_TOL:
            solutions.append((angle, result[1], result[2], result[0]))

    # 相邻有效采样偏差异号的区间二分求根
    for (a0, r0), (a1, r1) in zip(sampled, sampled[1:]):
        if r0 is None or r1 is None:
            continue  # 不跨无效采样连根
        d0, d1 = r0[0], r1[0]
        if d0 * d1 >= 0.0:
            continue  # 仅处理严格异号区间（零偏差采样已在扫描阶段记录）
        lo, hi = a0, a1
        d_lo = d0
        for _ in range(AIM_MAX_ITER):
            mid = 0.5 * (lo + hi)
            r_mid = evaluate(mid)
            if r_mid is None:
                break  # 二分遇无效追迹或无前向交点：丢弃该区间
            d_mid = r_mid[0]
            if abs(d_mid) <= AIM_DEV_TOL or (hi - lo) <= AIM_ANGLE_TOL_DEG:
                solutions.append((mid, r_mid[1], r_mid[2], d_mid))
                break
            if d_lo * d_mid < 0.0:
                hi = mid
            else:
                lo, d_lo = mid, d_mid

    # 按角度升序排序，按 1e-7 度去重（同簇保留 |偏差| 较小者）
    solutions.sort(key=lambda s: s[0])
    deduped: list[
        tuple[float, tuple[float, float], tuple[float, float], float]
    ] = []
    for sol in solutions:
        if deduped and sol[0] - deduped[-1][0] <= AIM_ANGLE_TOL_DEG:
            if abs(sol[3]) < abs(deduped[-1][3]):
                deduped[-1] = sol
            continue
        deduped.append(sol)
    return deduped


@dataclass
class _FocusPoint:
    """观察线沿法线平移到某个有符号位移处的光斑评估结果。

    ray_indices 为与平移线存在前向交点的光线编号（采样顺序，从 0 起）；
    命中少于两条时 centroid / rms 为 None，该位置无效。
    """

    shift: float
    hit_indices: tuple[int, ...]
    centroid: tuple[float, float] | None
    rms: float | None

    @property
    def valid(self) -> bool:
        return len(self.hit_indices) >= FOCUS_MIN_HITS


def focus_search(
    start: tuple[float, float],
    center_direction: tuple[float, float],
    half_angle_deg: float,
    sample_count: int,
    n_start: float,
    interfaces: list[tuple[tuple[float, float], tuple[float, float], float]],
    obs_point: tuple[float, float],
    obs_normal: tuple[float, float],
    shift_start: float,
    shift_end: float,
    scan_count: int,
) -> dict:
    """扇形采样光线只追迹一次，沿观察线单位法线等距平移观察线搜索 RMS 最低点。

    观察线平移：给定法线归一化为单位法线 n，位移 s 处的平移线经过
    obs_point + s·n（法线方向不变），s 从 shift_start 等间距取 scan_count
    个值到 shift_end（含端点）。每个采样位置统计**追迹完成**且与平移线有
    **前向**交点的光线，命中少于 FOCUS_MIN_HITS（2）条时该位置无效。

    对连续三个有效采样点，若中间点 RMS 不高于两端且三点命中光线编号集合
    相同，则视为一个低谷区间，用黄金分割法细化，最多 FOCUS_MAX_ROUNDS（24）
    轮（首轮取两个内点，之后每轮只新算一个内点）。细化试点一旦命中集合
    改变（或失效）即停止该区间，此前得到的有效试点全部保留。

    以所有有效扫描点和细化点中 RMS 最小者为最佳，平局取较小位移；没有
    任何有效位置时 best 为 None。该函数为即时计算，不写入数据库。

    返回 dict：
      ray_count / scan_count：扇形采样光线数与扫描点数；
      scan:   [(shift, hit_count, rms|None), ...]（按位移升序，含无效点）；
      refined:[(shift, hit_count, centroid, rms), ...]（有效细化点，按位移升序）；
      best:   (shift, centroid, rms) 或 None。
    """
    nrm = normalize(obs_normal)

    # 光线只追迹一次：缓存每条采样光线的出射段（未完成追迹的不纳入评估）
    ray_count = 0
    traced_rays: list[tuple[float, float]] = []  # 出射段 (origin, direction)
    for _angle, direction in sample_directions(
        center_direction, half_angle_deg, sample_count
    ):
        ray_count += 1
        segments, status, _reason = trace_ray(
            start, direction, n_start, interfaces
        )
        if status == "completed":
            traced_rays.append((segments[-1].origin, segments[-1].direction))

    cache: dict[float, _FocusPoint] = {}

    def evaluate(shift: float) -> _FocusPoint:
        """评估某个位移处的平移线光斑；同一位移只计算一次。"""
        if shift in cache:
            return cache[shift]
        moved = (obs_point[0] + shift * nrm[0], obs_point[1] + shift * nrm[1])
        hits: list[tuple[int, tuple[float, float]]] = []
        for idx, (origin, direction) in enumerate(traced_rays):
            point = _intersect_line(origin, direction, moved, nrm)
            if point is not None:
                hits.append((idx, point))
        indices = tuple(idx for idx, _p in hits)
        if len(hits) >= FOCUS_MIN_HITS:
            stats = centroid_and_rms([p for _idx, p in hits])
            centroid, rms = stats
        else:
            centroid = rms = None
        point = _FocusPoint(shift, indices, centroid, rms)
        cache[shift] = point
        return point

    # 等距扫描（含端点）
    scan_points: list[_FocusPoint] = []
    last = scan_count - 1
    for i in range(scan_count):
        shift = shift_start if i == 0 else shift_start + (shift_end - shift_start) * i / last
        scan_points.append(evaluate(shift))

    # 细化点池：(shift, _FocusPoint)，仅记录非扫描位置的有效试点；
    # 相邻低谷区间收缩时理论上可能探到同一位移，用集合保证每点只入池一次
    refined_points: list[tuple[float, _FocusPoint]] = []
    refined_shifts: set[float] = set()
    scan_shifts = {p.shift for p in scan_points}

    def probe(shift: float, base_set: frozenset[int]) -> _FocusPoint | None:
        """细化试点：返回有效且命中集合不变的点；失效或集合改变时返回 None。

        命中不少于两条的试点（无论集合是否改变）都作为有效点保留。
        """
        fp = evaluate(shift)
        if fp.valid and shift not in scan_shifts and shift not in refined_shifts:
            refined_shifts.add(shift)
            refined_points.append((shift, fp))
        if not fp.valid or frozenset(fp.hit_indices) != base_set:
            return None  # 命中集合改变（或失效）：停止该区间
        return fp

    # 连续三个有效扫描点构成低谷：rms[mid] <= rms[lo] 且 <= rms[hi]，命中集合相同
    for i in range(len(scan_points) - 2):
        p0, p1, p2 = scan_points[i], scan_points[i + 1], scan_points[i + 2]
        if not (p0.valid and p1.valid and p2.valid):
            continue  # 不跨无效扫描点
        if not (p1.rms <= p0.rms and p1.rms <= p2.rms):
            continue
        sets = frozenset(p0.hit_indices)
        if not (
            sets == frozenset(p1.hit_indices) and sets == frozenset(p2.hit_indices)
        ):
            continue

        lo, hi = p0.shift, p2.shift
        # 左内点 a = lo + 0.382L，右内点 b = lo + 0.618L
        a = lo + _GOLDEN_1MTAU * (hi - lo)
        b = lo + _GOLDEN_TAU * (hi - lo)
        fa = fb = None
        for round_i in range(FOCUS_MAX_ROUNDS):
            if round_i == 0:
                # 首轮取两个内点；任一试点命中集合改变即停止，不再算另一个
                fa = probe(a, sets)
                if fa is None:
                    break
                fb = probe(b, sets)
                if fb is None:
                    break
                continue
            # 后续每轮只新算一个内点，另一个由上轮幸存点复用
            if fa.rms <= fb.rms:
                # 最小点偏左：新区间 [lo, b]；幸存点 a 变为新左内点，新算右内点
                hi = b
                b, fb = a, fa
                a = lo + _GOLDEN_1MTAU * (hi - lo)
                fa = probe(a, sets)
            else:
                # 最小点偏右：新区间 [a, hi]；幸存点 b 变为新右内点，新算左内点
                lo = a
                a, fa = b, fb
                b = lo + _GOLDEN_TAU * (hi - lo)
                fb = probe(b, sets)
            if fa is None or fb is None:
                break

    # 以 RMS 最小为最佳，平局取较小位移
    candidates: list[tuple[float, _FocusPoint]] = [
        (p.shift, p) for p in scan_points if p.valid
    ]
    candidates.extend(refined_points)
    best: tuple[float, tuple[float, float], float] | None = None
    if candidates:
        best_shift, best_fp = min(
            candidates, key=lambda item: (item[1].rms, item[0])
        )
        best = (best_shift, best_fp.centroid, best_fp.rms)

    refined_points.sort(key=lambda item: item[0])
    return {
        "ray_count": ray_count,
        "scan_count": scan_count,
        "scan": [(p.shift, len(p.hit_indices), p.rms) for p in scan_points],
        "refined": [
            (p.shift, len(p.hit_indices), p.centroid, p.rms)
            for _shift, p in refined_points
        ],
        "best": best,
    }
