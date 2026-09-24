"""二维介质光线追迹 API。

POST /api/traces                      提交光线与界面序列，返回追迹结果并保存
GET  /api/traces                      列出所有已保存的追迹（摘要）
GET  /api/traces/{trace_id}           按编号查询完整追迹（含各段交点与方向）
GET  /api/traces/{trace_id}/segments/{seg_index}  按编号查询某一段
POST /api/beam-evaluations            以中心方向为轴对称采样多条光线，返回每条光线
                                      的发射角、追迹状态、观察线交点及质心/RMS 半径
POST /api/power-evaluations           复用光束采样与追迹，按菲涅耳公式逐界面累乘
                                      s/p 透过率，返回各光线功率与功率加权质心
POST /api/aim-solutions               在半角范围内扫描发射角并对偏差异号区间二分，
                                      求解命中观察线目标偏移量的光线，返回角度、
                                      出射方向、观察线交点与目标偏差（不入库）
POST /api/alignment-sensitivities     基准追迹后逐界面将法线旋转 ±扰动角重新追迹，
                                      按中心差分求观察线偏移对法线角度的导数（不入库）
"""
from __future__ import annotations

from contextlib import asynccontextmanager

from fastapi import FastAPI, HTTPException

from . import db
from .geometry import (
    aim_solutions,
    alignment_sensitivities,
    centroid_and_rms,
    observation_intersection,
    ray_transmittance,
    sample_directions,
    trace_ray,
    weighted_centroid,
)
from .schemas import (
    AimRequest,
    AimSolutionOut,
    AlignmentSensitivityOut,
    BeamEvaluationOut,
    BeamRequest,
    PowerEvaluationOut,
    PowerRequest,
    SegmentOut,
    SensitivityRequest,
    TraceOut,
    TraceRequest,
    TraceSummary,
)


@asynccontextmanager
async def lifespan(app: FastAPI):
    db.init_db()  # 启动时建库建表
    yield


app = FastAPI(
    title="2D Medium Ray Tracing API",
    version="1.0.0",
    description="二维介质光线追迹：按斯涅尔定律折射，遇全反射停止并返回原因。",
    lifespan=lifespan,
)


@app.get("/")
def root() -> dict:
    return {"service": "2d-ray-tracing", "docs": "/docs", "api": "/api/traces"}


@app.post("/api/traces", response_model=TraceOut, status_code=201)
def create_trace(req: TraceRequest) -> dict:
    interfaces = [
        (itf.point.as_tuple(), itf.normal.as_tuple(), itf.n_next)
        for itf in req.interfaces
    ]
    segments, status, reason = trace_ray(
        req.start.as_tuple(), req.direction.as_tuple(), req.n_start, interfaces
    )
    trace_id = db.insert_trace(
        req.start.as_tuple(),
        req.direction.as_tuple(),
        req.n_start,
        status,
        reason,
        segments,
    )
    return db.get_trace(trace_id)


@app.get("/api/traces", response_model=list[TraceSummary])
def list_traces() -> list[dict]:
    return db.list_traces()


@app.get("/api/traces/{trace_id}", response_model=TraceOut)
def get_trace(trace_id: int) -> dict:
    trace = db.get_trace(trace_id)
    if trace is None:
        raise HTTPException(status_code=404, detail=f"追迹 {trace_id} 不存在")
    return trace


@app.get("/api/traces/{trace_id}/segments/{seg_index}", response_model=SegmentOut)
def get_segment(trace_id: int, seg_index: int) -> dict:
    segment = db.get_segment(trace_id, seg_index)
    if segment is None:
        raise HTTPException(
            status_code=404, detail=f"追迹 {trace_id} 的第 {seg_index} 段不存在"
        )
    return segment


@app.post("/api/beam-evaluations", response_model=BeamEvaluationOut)
def evaluate_beam(req: BeamRequest) -> dict:
    """以中心方向为 0° 在 [-a, +a] 内等间距旋转采样，逐条追迹并评估观察线光斑。

    仅对追迹完成（status=completed）且与观察线存在前向交点的光线纳入统计；
    未完成、与观察线平行或交点在光线反方向时 observation_point 为 null。
    """
    interfaces = [
        (itf.point.as_tuple(), itf.normal.as_tuple(), itf.n_next)
        for itf in req.interfaces
    ]
    start = req.start.as_tuple()

    rays_out: list[dict] = []
    hit_points: list[tuple[float, float]] = []
    for angle, direction in sample_directions(
        req.center_direction.as_tuple(), req.half_angle_deg, req.sample_count
    ):
        segments, status, reason = trace_ray(
            start, direction, req.n_start, interfaces
        )
        observation_point = None
        if status == "completed":
            observation_point = observation_intersection(
                segments,
                req.observation_line.point.as_tuple(),
                req.observation_line.normal.as_tuple(),
            )
            if observation_point is not None:
                hit_points.append(observation_point)
        rays_out.append(
            {
                "emission_angle_deg": angle,
                "direction": {"x": direction[0], "y": direction[1]},
                "status": status,
                "reason": reason,
                "observation_point": (
                    None
                    if observation_point is None
                    else {"x": observation_point[0], "y": observation_point[1]}
                ),
            }
        )

    stats = centroid_and_rms(hit_points)
    return {
        "half_angle_deg": req.half_angle_deg,
        "sample_count": req.sample_count,
        "ray_count": len(rays_out),
        "hit_count": len(hit_points),
        "centroid": None if stats is None else {"x": stats[0][0], "y": stats[0][1]},
        "rms_radius": None if stats is None else stats[1],
        "rays": rays_out,
    }


@app.post("/api/power-evaluations", response_model=PowerEvaluationOut)
def evaluate_power(req: PowerRequest) -> dict:
    """复用光束评估的扇形采样与光线追迹，按菲涅耳公式评估 s/p 偏振透过功率。

    每条光线入射功率为 1，s/p 偏振分别占 s_fraction 与 1 - s_fraction；
    逐界面由方向与法线求 cos_i、cos_t，按
        Rs=[(n1 cos_i-n2 cos_t)/(n1 cos_i+n2 cos_t)]²、
        Rp=[(n2 cos_i-n1 cos_t)/(n2 cos_i+n1 cos_t)]²
    累乘透过率 1 - R。全反射或漏过界面（追迹未完成）记零功率；仅追迹完成
    且与观察线有前向交点的光线计入有效功率与功率加权质心，有效功率为零时
    质心为 null。该接口为即时计算，结果不写入数据库。
    """
    interfaces = [
        (itf.point.as_tuple(), itf.normal.as_tuple(), itf.n_next)
        for itf in req.interfaces
    ]
    start = req.start.as_tuple()

    rays_out: list[dict] = []
    weighted_hits: list[tuple[tuple[float, float], float]] = []
    effective_power = 0.0
    for angle, direction in sample_directions(
        req.center_direction.as_tuple(), req.half_angle_deg, req.sample_count
    ):
        segments, status, reason = trace_ray(
            start, direction, req.n_start, interfaces
        )
        t_s, t_p = ray_transmittance(segments, status, interfaces)
        power = req.s_fraction * t_s + (1.0 - req.s_fraction) * t_p

        observation_point = None
        if status == "completed":
            observation_point = observation_intersection(
                segments,
                req.observation_line.point.as_tuple(),
                req.observation_line.normal.as_tuple(),
            )
            if observation_point is not None:
                effective_power += power
                weighted_hits.append((observation_point, power))

        rays_out.append(
            {
                "emission_angle_deg": angle,
                "direction": {"x": direction[0], "y": direction[1]},
                "status": status,
                "reason": reason,
                "transmittance_s": t_s,
                "transmittance_p": t_p,
                "power": power,
                "observation_point": (
                    None
                    if observation_point is None
                    else {"x": observation_point[0], "y": observation_point[1]}
                ),
            }
        )

    centroid = weighted_centroid(weighted_hits)
    return {
        "half_angle_deg": req.half_angle_deg,
        "sample_count": req.sample_count,
        "s_fraction": req.s_fraction,
        "ray_count": len(rays_out),
        "hit_count": len(weighted_hits),
        "total_input_power": float(len(rays_out)),
        "effective_power": effective_power,
        "weighted_centroid": (
            None if centroid is None else {"x": centroid[0], "y": centroid[1]}
        ),
        "rays": rays_out,
    }


@app.post("/api/aim-solutions", response_model=list[AimSolutionOut])
def compute_aim_solutions(req: AimRequest) -> list[dict]:
    """扫描 + 二分求解使出射光线命中观察线目标位置的发射角。

    目标偏移量：沿观察线单位切向量（法线逆时针旋转 90°）相对观察线给定点的
    有符号距离。在半角范围内等间距扫描，仅对相邻有效采样偏差异号的区间二分；
    采样点直接命中（|偏差| ≤ 1e-6）也记为解。结果按角度升序、按 1e-7 度去重，
    无解时返回空列表。该接口为即时计算，不写入数据库。
    """
    interfaces = [
        (itf.point.as_tuple(), itf.normal.as_tuple(), itf.n_next)
        for itf in req.interfaces
    ]
    solutions = aim_solutions(
        req.start.as_tuple(),
        req.center_direction.as_tuple(),
        req.half_angle_deg,
        req.sample_count,
        req.n_start,
        interfaces,
        req.observation_line.point.as_tuple(),
        req.observation_line.normal.as_tuple(),
        req.target_offset,
    )
    return [
        {
            "emission_angle_deg": angle,
            "exit_direction": {"x": exit_dir[0], "y": exit_dir[1]},
            "observation_point": {"x": point[0], "y": point[1]},
            "target_deviation": deviation,
        }
        for angle, exit_dir, point, deviation in solutions
    ]


def _sensitivity_run(run: tuple, normal: tuple[float, float] | None = None,
                     angle: float | None = None) -> dict:
    """将几何层一次追迹结果 (status, reason, point, offset[, nrm]) 转为响应字典。"""
    status, reason, point, offset = run[:4]
    out = {
        "status": status,
        "reason": reason,
        "observation_point": (
            None if point is None else {"x": point[0], "y": point[1]}
        ),
        "offset": offset,
    }
    if normal is not None:
        out["perturbed_normal"] = {"x": normal[0], "y": normal[1]}
        out["perturbation_angle_deg"] = angle
    return out


@app.post("/api/alignment-sensitivities", response_model=AlignmentSensitivityOut)
def compute_alignment_sensitivities(req: SensitivityRequest) -> dict:
    """基准追迹后逐界面扰动法线 ±δ，按中心差分求观察线偏移对法线角度的导数。

    观察线法线归一化后逆时针旋转 90° 得单位切向量；先计算基准出射线前向交点
    相对观察线给定点的有符号偏移。对每个界面只旋转其法线 ±perturbation_angle_deg
    （界面点及其他参数不变）重新追迹，导数 = (偏移(+δ) - 偏移(-δ)) / (2δ)，
    单位为长度/度。基准或任一扰动发生全反射、漏过界面或无法前向命中观察线时
    （status 非 completed），该次偏移为 null、对应导数为 null，不跨失效分支
    求差分。该接口为即时计算，不写入数据库。
    """
    interfaces = [
        (itf.point.as_tuple(), itf.normal.as_tuple(), itf.n_next)
        for itf in req.interfaces
    ]
    result = alignment_sensitivities(
        req.start.as_tuple(),
        req.direction.as_tuple(),
        req.n_start,
        interfaces,
        req.observation_line.point.as_tuple(),
        req.observation_line.normal.as_tuple(),
        req.perturbation_angle_deg,
    )
    tangent = result["tangent"]
    interfaces_out: list[dict] = []
    for idx, (r_plus, r_minus, derivative) in enumerate(result["interfaces"]):
        interfaces_out.append(
            {
                "interface_index": idx,
                "plus": _sensitivity_run(
                    r_plus, r_plus[4], req.perturbation_angle_deg
                ),
                "minus": _sensitivity_run(
                    r_minus, r_minus[4], -req.perturbation_angle_deg
                ),
                "derivative_per_deg": derivative,
            }
        )
    return {
        "tangent": {"x": tangent[0], "y": tangent[1]},
        "perturbation_angle_deg": req.perturbation_angle_deg,
        "baseline": _sensitivity_run(result["baseline"]),
        "interfaces": interfaces_out,
    }
