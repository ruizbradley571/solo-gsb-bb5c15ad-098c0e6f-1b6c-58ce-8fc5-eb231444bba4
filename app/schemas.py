"""请求 / 响应数据模型。"""
from __future__ import annotations

import math

from pydantic import BaseModel, Field, field_validator

_EPS = 1e-12


class Vec2(BaseModel):
    x: float
    y: float

    def as_tuple(self) -> tuple[float, float]:
        return (self.x, self.y)


class InterfaceIn(BaseModel):
    """一个界面：线上一点、指向下一介质的法线、下一介质折射率。"""

    point: Vec2
    normal: Vec2
    n_next: float = Field(gt=0, description="下一介质折射率，必须为正")

    @field_validator("normal")
    @classmethod
    def normal_nonzero(cls, v: Vec2) -> Vec2:
        if math.hypot(v.x, v.y) < _EPS:
            raise ValueError("法线不能为零向量")
        return v


class TraceRequest(BaseModel):
    start: Vec2
    direction: Vec2
    n_start: float = Field(gt=0, description="起始介质折射率，必须为正")
    interfaces: list[InterfaceIn] = Field(default_factory=list)

    @field_validator("direction")
    @classmethod
    def direction_nonzero(cls, v: Vec2) -> Vec2:
        if math.hypot(v.x, v.y) < _EPS:
            raise ValueError("方向不能为零向量")
        return v


class SegmentOut(BaseModel):
    seg_index: int
    origin: Vec2
    direction: Vec2
    n_medium: float
    hit_point: Vec2 | None
    event: str | None
    n_after: float | None


class TraceOut(BaseModel):
    id: int
    start: Vec2
    direction: Vec2
    n_start: float
    status: str
    reason: str | None
    created_at: str
    segments: list[SegmentOut]


class TraceSummary(BaseModel):
    id: int
    status: str
    reason: str | None
    n_start: float
    created_at: str


class ObservationLineIn(BaseModel):
    """观察线：线上一点与法线（法线方向不影响直线位置，仅用于求交方程）。"""

    point: Vec2
    normal: Vec2

    @field_validator("normal")
    @classmethod
    def normal_nonzero(cls, v: Vec2) -> Vec2:
        if math.hypot(v.x, v.y) < _EPS:
            raise ValueError("法线不能为零向量")
        return v


class BeamRequest(BaseModel):
    start: Vec2
    center_direction: Vec2
    half_angle_deg: float = Field(gt=0, lt=90, description="对称半角（度），0<a<90")
    sample_count: int = Field(ge=2, le=101, description="采样数，2～101（含端点）")
    n_start: float = Field(gt=0, description="起始介质折射率，必须为正")
    interfaces: list[InterfaceIn] = Field(default_factory=list)
    observation_line: ObservationLineIn

    @field_validator("center_direction")
    @classmethod
    def direction_nonzero(cls, v: Vec2) -> Vec2:
        if math.hypot(v.x, v.y) < _EPS:
            raise ValueError("中心方向不能为零向量")
        return v


class BeamRayOut(BaseModel):
    emission_angle_deg: float
    direction: Vec2
    status: str
    reason: str | None
    observation_point: Vec2 | None


class BeamEvaluationOut(BaseModel):
    half_angle_deg: float
    sample_count: int
    ray_count: int
    hit_count: int
    centroid: Vec2 | None
    rms_radius: float | None
    rays: list[BeamRayOut]


class PowerRequest(BeamRequest):
    """功率评估请求：在光束请求基础上增加 s 偏振入射功率占比。

    每条光线入射功率为 1，其中 s 偏振占 s_fraction、p 偏振占 1 - s_fraction。
    """

    s_fraction: float = Field(ge=0, le=1, description="s 偏振入射功率占比，取值 0～1")


class PowerRayOut(BaseModel):
    emission_angle_deg: float
    direction: Vec2
    status: str
    reason: str | None
    transmittance_s: float
    transmittance_p: float
    power: float
    observation_point: Vec2 | None


class PowerEvaluationOut(BaseModel):
    half_angle_deg: float
    sample_count: int
    s_fraction: float
    ray_count: int
    hit_count: int
    total_input_power: float
    effective_power: float
    weighted_centroid: Vec2 | None
    rays: list[PowerRayOut]


class AimRequest(BeamRequest):
    """瞄准求解请求：在光束请求基础上增加目标偏移量。

    目标偏移量：沿观察线单位切向量（法线逆时针旋转 90°）相对观察线给定点的
    有符号距离，即目标点 = 观察线给定点 + target_offset · 切向量。
    """

    target_offset: float = Field(
        description="目标点沿观察线单位切向量相对观察线给定点的有符号距离"
    )


class AimSolutionOut(BaseModel):
    emission_angle_deg: float
    exit_direction: Vec2
    observation_point: Vec2
    target_deviation: float


class SensitivityRequest(BaseModel):
    """法线灵敏度请求：单条光线、顺序界面、观察线与正值扰动角（度）。"""

    start: Vec2
    direction: Vec2
    n_start: float = Field(gt=0, description="起始介质折射率，必须为正")
    interfaces: list[InterfaceIn] = Field(default_factory=list)
    observation_line: ObservationLineIn
    perturbation_angle_deg: float = Field(
        gt=0, description="法线扰动角（度），必须为正；对每个界面分别旋转 ±该角度"
    )

    @field_validator("direction")
    @classmethod
    def direction_nonzero(cls, v: Vec2) -> Vec2:
        if math.hypot(v.x, v.y) < _EPS:
            raise ValueError("方向不能为零向量")
        return v


class SensitivityRunOut(BaseModel):
    """一次追迹（基准或单侧扰动）的状态与观察线偏移结果。"""

    status: str
    reason: str | None
    observation_point: Vec2 | None
    offset: float | None


class SensitivityPerturbationOut(SensitivityRunOut):
    """单侧法线扰动结果：附本次扰动后的单位法线与扰动角（度，逆时针为正）。"""

    perturbed_normal: Vec2
    perturbation_angle_deg: float


class InterfaceSensitivityOut(BaseModel):
    """单个界面的两次法线扰动结果与中心差分导数（长度/度）。"""

    interface_index: int
    plus: SensitivityPerturbationOut
    minus: SensitivityPerturbationOut
    derivative_per_deg: float | None


class AlignmentSensitivityOut(BaseModel):
    tangent: Vec2
    perturbation_angle_deg: float
    baseline: SensitivityRunOut
    interfaces: list[InterfaceSensitivityOut]
