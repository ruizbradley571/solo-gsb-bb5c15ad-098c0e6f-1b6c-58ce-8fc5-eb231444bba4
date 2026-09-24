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

FOCUS_GOLDEN_MAX_ROUNDS = 24  # 单个低谷区间黄金分割细化的最多轮数
GOLDEN_RATIO = (math.sqrt(5.0) - 1.0) / 2.0  # (√5-1)/2 ≈ 0.618


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


def _observation_point_on_line(
    origin: tuple[float, float],
    direction: tuple[float, float],
    line_point: tuple[float, float],
    line_normal: tuple[float, float],
) -> tuple[float, float] | None:
    """给定出射段原点/方向，求其与指定直线的前向交点；平行或交点在身后返回 None。"""
    return _intersect_line(origin, direction, line_point, line_normal)


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
    scan_points: int,
) -> dict:
    """沿观察线单位法线平移观察线做焦面扫描，并用黄金分割细化 RMS 低谷。

    光线只追迹一次：先按光束评估相同的扇形采样逐条追迹，保留各光线的出射段；
    随后将观察线（点 obs_point、单位法线 n）沿 n 平移位移 s（直线点变为
    obs_point + s·n，法线不变），逐扫描位置 s 求出射段与平移线的前向交点。

    - 仅追迹完成（status=completed）且与平移线有前向交点的光线命中；命中少于
      两条（命中编号集合大小 < 2）时该位置无效（hit_count/centroid/rms 为 None）。
    - 在位移 [shift_start, shift_end]（shift_start < shift_end）内含端点等间距取
      scan_points 个扫描点。
    - 对连续三个有效扫描点，中间点 RMS 不高于两端（r1 ≤ r0 且 r1 ≤ r2）且三点
      命中光线编号集合相同的低谷区间，以黄金分割法细化，最多 24 轮；细化中一旦
      试点的命中集合改变（与区间端点集合不同）或试点无效，立即停止该区间，已得
      到的有效试点保留。
    - 以所有有效扫描点与细化点中 RMS 最小者为最佳，平局取较小位移；无任何有效
      位置时 best 为 None。

    返回 dict：
      scan:     [(shift, hit_indices|None, centroid|None, rms|None), ...]
      refined:  [(shift, hit_indices, centroid, rms), ...]（有效细化试点）
      best:     (shift, centroid, rms) | None
    """
    nrm = normalize(obs_normal)

    # 光线只追迹一次：记录每条完成追迹光线的出射段（原点、方向）与其编号
    fan: list[tuple[int, tuple[float, float], tuple[float, float]]] = []
    for idx, (_angle, direction) in enumerate(
        sample_directions(center_direction, half_angle_deg, sample_count)
    ):
        segments, status, _reason = trace_ray(start, direction, n_start, interfaces)
        if status == "completed":
            fan.append((idx, segments[-1].origin, segments[-1].direction))

    def probe(
        shift: float,
    ) -> tuple[frozenset[int], tuple[float, float], float] | None:
        """在沿法线平移 shift 的观察线上求命中编号集合、质心与 RMS；无效返回 None。"""
        line_point = (
            obs_point[0] + shift * nrm[0],
            obs_point[1] + shift * nrm[1],
        )
        hits: list[tuple[int, tuple[float, float]]] = []
        for ray_idx, seg_origin, seg_dir in fan:
            point = _observation_point_on_line(seg_origin, seg_dir, line_point, nrm)
            if point is not None:
                hits.append((ray_idx, point))
        if len(hits) < 2:
            return None  # 少于两条命中：该位置无效
        points = [point for _idx, point in hits]
        centroid, rms = centroid_and_rms(points)
        return frozenset(idx for idx, _p in hits), centroid, rms

    # 等间距扫描（端点精确取 shift_start / shift_end）
    shifts = [
        shift_start
        if i == 0
        else shift_end
        if i == scan_points - 1
        else shift_start
        + (shift_end - shift_start) * i / (scan_points - 1)
        for i in range(scan_points)
    ]
    scan: list[
        tuple[float, frozenset[int] | None, tuple[float, float] | None, float | None]
    ] = []
    evaluated: list[
        tuple[float, frozenset[int], tuple[float, float], float]
    ] = []
    for shift in shifts:
        result = probe(shift)
        if result is None:
            scan.append((shift, None, None, None))
        else:
            hit_set, centroid, rms = result
            scan.append((shift, hit_set, centroid, rms))
            evaluated.append((shift, hit_set, centroid, rms))

    refined: list[
        tuple[float, frozenset[int], tuple[float, float], float]
    ] = []

    def valid_scan(i: int) -> bool:
        return scan[i][2] is not None

    # 连续三个有效扫描点构成的低谷区间逐一黄金分割细化
    for i in range(scan_points - 2):
        if not (valid_scan(i) and valid_scan(i + 1) and valid_scan(i + 2)):
            continue
        s0, set0, _c0, r0 = scan[i]
        s1, set1, _c1, r1 = scan[i + 1]
        s2, set2, _c2, r2 = scan[i + 2]
        if not (r1 <= r0 and r1 <= r2):
            continue  # 中间点不高于两端才构成低谷
        if not (set0 == set1 == set2):
            continue  # 命中光线编号必须相同
        target_set = set0
        lo, hi = float(s0), float(s2)

        def trial(
            s: float,
        ) -> tuple[frozenset[int], tuple[float, float], float] | None:
            """细化试点：命中集合改变或无效时返回 None（调用方据此停止区间）。"""
            result = probe(s)
            if result is None:
                return None
            hit_set, centroid, rms = result
            if hit_set != target_set:
                return None  # 命中集合改变：停止该区间
            return result

        width = hi - lo
        c = lo + (1.0 - GOLDEN_RATIO) * width
        d = lo + GOLDEN_RATIO * width
        rc = trial(c)
        if rc is None:
            continue  # 首个内部试点即失效：区间无有效细化点
        hit_c, cent_c, r_c = rc
        refined.append((c, hit_c, cent_c, r_c))
        evaluated.append((c, hit_c, cent_c, r_c))
        rd = trial(d)
        if rd is None:
            continue  # 保留 c 处已得到的有效试点，停止该区间
        hit_d, cent_d, r_d = rd
        refined.append((d, hit_d, cent_d, r_d))
        evaluated.append((d, hit_d, cent_d, r_d))

        for _round in range(FOCUS_GOLDEN_MAX_ROUNDS):
            if r_c < r_d:
                hi, d, r_d, hit_d, cent_d = d, c, r_c, hit_c, cent_c
                width = hi - lo
                c = lo + (1.0 - GOLDEN_RATIO) * width
                rc_new = trial(c)
                if rc_new is None:
                    break  # 命中集合改变或试点无效：停止该区间，保留已得试点
                hit_c, cent_c, r_c = rc_new
                refined.append((c, hit_c, cent_c, r_c))
                evaluated.append((c, hit_c, cent_c, r_c))
            else:
                lo, c, r_c, hit_c, cent_c = c, d, r_d, hit_d, cent_d
                width = hi - lo
                d = lo + GOLDEN_RATIO * width
                rd_new = trial(d)
                if rd_new is None:
                    break
                hit_d, cent_d, r_d = rd_new
                refined.append((d, hit_d, cent_d, r_d))
                evaluated.append((d, hit_d, cent_d, r_d))

    # RMS 最小者为最佳，平局取较小位移
    best: tuple[float, tuple[float, float], float] | None = None
    for shift, _hit_set, centroid, rms in evaluated:
        if best is None or rms < best[2] or (rms == best[2] and shift < best[0]):
            best = (shift, centroid, rms)

    return {"scan": scan, "refined": refined, "best": best}
