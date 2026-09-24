# 2D 介质光线追迹 API

基于 **Python + FastAPI + SQLite** 的二维光线追迹服务：提交初始光线（起点、方向）、起始介质折射率，以及按顺序排列的界面（界面上一点、指向下一介质的法线、下一介质折射率），服务逐界面计算前向交点并依斯涅尔定律更新方向；遇到全反射（或光线与界面平行 / 交点在反方向）即停止并返回原因。追迹结果保存到 SQLite，可按编号查询各段交点与方向。

## 启动

```bash
docker compose up --build
```

启动时会安装依赖并自动建库建表（SQLite 文件存于命名卷 `raytrace-data` 中）。

- **API 地址**：http://localhost:8000
- **交互式文档（Swagger UI）**：http://localhost:8000/docs
- **端口**：`8000`

本地开发（不用 Docker）：

```bash
pip install -r requirements.txt
uvicorn app.main:app --host 0.0.0.0 --port 8000
```

## API 一览

| 方法 | 路径 | 说明 |
| --- | --- | --- |
| POST | `/api/traces` | 提交光线与界面序列，执行追迹、保存并返回完整结果 |
| GET | `/api/traces` | 列出所有已保存追迹的摘要 |
| GET | `/api/traces/{trace_id}` | 按编号查询完整追迹（含各段交点与方向） |
| GET | `/api/traces/{trace_id}/segments/{seg_index}` | 按编号查询某一段的交点与方向 |
| POST | `/api/beam-evaluations` | 扇形采样多条光线逐条追迹，返回各光线发射角、状态、观察线交点及质心/RMS（不入库） |
| POST | `/api/power-evaluations` | 复用光束采样与追迹，按菲涅耳公式逐界面累乘 s/p 透过率，返回各光线功率与功率加权质心（不入库） |
| POST | `/api/aim-solutions` | 扫描发射角并对偏差异号区间二分，求解命中观察线目标偏移量的光线（不入库） |
| POST | `/api/alignment-sensitivities` | 基准追迹后逐界面将法线旋转 ±扰动角重新追迹，中心差分求观察线偏移对法线角度的导数（不入库） |

### 请求示例

```bash
curl -X POST http://localhost:8000/api/traces \
  -H "Content-Type: application/json" \
  -d '{
    "start": {"x": 0, "y": 0},
    "direction": {"x": 1, "y": 0.5},
    "n_start": 1.0,
    "interfaces": [
      {"point": {"x": 2, "y": 0}, "normal": {"x": 1, "y": 0}, "n_next": 1.5},
      {"point": {"x": 4, "y": 0}, "normal": {"x": 1, "y": 0}, "n_next": 1.0}
    ]
  }'
```

### 响应说明

- `status`：`completed`（穿过所有界面）、`total_internal_reflection`（全反射停止）、`no_intersection`（与界面平行或交点在反方向）
- `reason`：停止原因（完成时为 `null`）
- `segments[]`：每段包含 `origin`（段起点）、`direction`（段方向）、`n_medium`（所在介质折射率）、`hit_point`（与下一界面的交点，出射段为 `null`）、`event`（`refraction` / `total_internal_reflection` / `no_intersection` / `null`）

### 全反射示例

```bash
curl -X POST http://localhost:8000/api/traces \
  -H "Content-Type: application/json" \
  -d '{
    "start": {"x": 0, "y": 0},
    "direction": {"x": 1, "y": 2},
    "n_start": 1.5,
    "interfaces": [
      {"point": {"x": 2, "y": 0}, "normal": {"x": 1, "y": 0}, "n_next": 1.0}
    ]
  }'
# status = "total_internal_reflection"，reason 说明入射角超过临界角
```

## 光束评估

`POST /api/beam-evaluations` 以中心方向为 0°，在对称半角 `[-half_angle_deg, +half_angle_deg]`（度，要求 `0 < a < 90`）内**含端点等间距**旋转生成 `sample_count`（2～101）条光线，逐条复用同一套追迹逻辑穿过给定的顺序界面，并在追迹完成后取最终出射段与**观察线**求前向交点。该接口为即时计算，结果不写入数据库。

- 发射角正方向为数学上的**逆时针**方向，端点光线的发射角精确为 `-a` 与 `+a`。
- 仅当追迹状态为 `completed` 且出射段与观察线存在**前向**交点时返回交点；追迹未完成（全反射 / 中途无交点）、与观察线平行或交点位于光线反方向时 `observation_point` 为 `null`。
- `centroid` / `rms_radius` 只统计非空交点：质心为交点坐标均值，RMS 半径为各交点到质心距离平方均值的开方；没有任何命中时二者均为 `null`。

| 方法 | 路径 | 说明 |
| --- | --- | --- |
| POST | `/api/beam-evaluations` | 光束扇形采样、逐条追迹并评估观察线上的光斑（质心与 RMS 半径） |

### 光束评估请求示例

```bash
curl -X POST http://localhost:8000/api/beam-evaluations \
  -H "Content-Type: application/json" \
  -d '{
    "start": {"x": 0, "y": 0},
    "center_direction": {"x": 1, "y": 0},
    "half_angle_deg": 10,
    "sample_count": 21,
    "n_start": 1.0,
    "interfaces": [
      {"point": {"x": 2, "y": 0}, "normal": {"x": 1, "y": 0}, "n_next": 1.5},
      {"point": {"x": 4, "y": 0}, "normal": {"x": 1, "y": 0}, "n_next": 1.0}
    ],
    "observation_line": {"point": {"x": 10, "y": 0}, "normal": {"x": 1, "y": 0}}
  }'
```

### 光束评估响应说明

```json
{
  "half_angle_deg": 10.0,
  "sample_count": 21,
  "ray_count": 21,
  "hit_count": 21,
  "centroid": {"x": 10.0, "y": 0.0},
  "rms_radius": 0.9922,
  "rays": [
    {
      "emission_angle_deg": -10.0,
      "direction": {"x": 0.9848, "y": -0.1736},
      "status": "completed",
      "reason": null,
      "observation_point": {"x": 10.0, "y": -0.9922}
    }
  ]
}
```

- `ray_count`：实际追迹的光线数（等于 `sample_count`）；`hit_count`：观察线上存在有效前向交点的光线数。
- `rays[]` 按发射角从 `-a` 到 `+a` 排列，每条含 `emission_angle_deg`（相对中心方向的发射角，度）、`direction`（归一化后的实际发射方向）、`status` / `reason`（复用单光线追迹的状态与停止原因）、`observation_point`（观察线前向交点，无效时为 `null`）。
- 无任何命中时：`hit_count` 为 0，`centroid` 与 `rms_radius` 均为 `null`。

## 功率评估

`POST /api/power-evaluations` 在光束请求（起点、中心方向、折射率、顺序界面、观察线、半角、采样数）基础上增加 `s_fraction`（s 偏振入射功率占比，取值 0～1，p 偏振占 `1 - s_fraction`），复用同一套扇形采样与光线追迹，逐界面按**菲涅耳公式**计算 s/p 偏振的功率反射率与透过率。该接口为即时计算，结果不写入数据库。

- 以光线方向与界面法线求入射角余弦 `cos_i`（光线从法线背面入射时翻转法线），由斯涅尔定律 `sin_t = (n1/n2) sin_i` 求折射角余弦 `cos_t`：
  - `Rs = [(n1 cos_i - n2 cos_t)/(n1 cos_i + n2 cos_t)]²`
  - `Rp = [(n2 cos_i - n1 cos_t)/(n2 cos_i + n1 cos_t)]²`
  - 逐界面累乘**功率透过率** `1 - R`：`Ts = ∏(1 - Rs)`、`Tp = ∏(1 - Rp)`。
- 每条光线入射功率为 **1**，按占比混合：`power = s_fraction · Ts + (1 - s_fraction) · Tp`。
- **全反射**（任一束光线追迹状态为 `total_internal_reflection`）或**漏过界面**（`no_intersection`，与界面平行或交点在反方向）时该光线记**零功率**（`Ts = Tp = power = 0`）。
- 仅**追迹完成且与观察线有前向交点**的光线计入**有效功率**（其 `power` 之和）与**功率加权质心** `Σ(power·point)/Σpower`；追迹完成但观察线无交点的光线有透过率但不计入。
- `total_input_power` 为总入射功率（等于光线数，每条入射功率为 1）；有效功率为零时 `weighted_centroid` 为 `null`。

| 方法 | 路径 | 说明 |
| --- | --- | --- |
| POST | `/api/power-evaluations` | 扇形采样、逐条追迹并计算 s/p 菲涅耳透过率、有效功率与功率加权质心 |

### 功率评估请求示例

```bash
curl -X POST http://localhost:8000/api/power-evaluations \
  -H "Content-Type: application/json" \
  -d '{
    "start": {"x": 0, "y": 0},
    "center_direction": {"x": 1, "y": 0},
    "half_angle_deg": 10,
    "sample_count": 5,
    "n_start": 1.0,
    "interfaces": [
      {"point": {"x": 2, "y": 0}, "normal": {"x": 1, "y": 0}, "n_next": 1.5},
      {"point": {"x": 4, "y": 0}, "normal": {"x": 1, "y": 0}, "n_next": 1.0}
    ],
    "observation_line": {"point": {"x": 10, "y": 0}, "normal": {"x": 1, "y": 0}},
    "s_fraction": 0.5
  }'
```

### 功率评估响应说明

```json
{
  "half_angle_deg": 10.0,
  "sample_count": 5,
  "s_fraction": 0.5,
  "ray_count": 5,
  "hit_count": 5,
  "total_input_power": 5.0,
  "effective_power": 4.607942624633463,
  "weighted_centroid": {"x": 10.0, "y": 0.0},
  "rays": [
    {
      "emission_angle_deg": -10.0,
      "direction": {"x": 0.984807753012208, "y": -0.17364817766693033},
      "status": "completed",
      "reason": null,
      "transmittance_s": 0.918416539493318,
      "transmittance_p": 0.9247294040133414,
      "power": 0.9215729717533296,
      "observation_point": {"x": 10.0, "y": -1.6437139662827691}
    },
    {
      "emission_angle_deg": 0.0,
      "direction": {"x": 1.0, "y": 0.0},
      "status": "completed",
      "reason": null,
      "transmittance_s": 0.9216,
      "transmittance_p": 0.9216,
      "power": 0.9216,
      "observation_point": {"x": 10.0, "y": 0.0}
    }
  ]
}
```

- `total_input_power`：总入射功率，恒等于 `ray_count`（每条光线入射功率为 1）。
- `effective_power`：有效功率，仅累加追迹完成且与观察线有前向交点的光线的 `power`。
- `weighted_centroid`：有效光线按 `power` 加权的质心；有效功率为 0 时为 `null`。
- `rays[].transmittance_s` / `transmittance_p`：该光线逐界面累乘的 s/p 总透过率（`0～1`）；全反射或漏过界面时为 `0`；无界面时为 `1`。
- `rays[].power`：按 `s_fraction` 混合后的该光线功率；`observation_point` 为观察线前向交点，无效时为 `null`。

## 瞄准求解

`POST /api/aim-solutions` 在光束请求（起点、中心方向、折射率、顺序界面、观察线、半角、采样数）的基础上增加 `target_offset`（目标偏移量），求解使出射光线命中观察线上目标位置的发射角。该接口为即时计算，结果不写入数据库。

- **目标偏移量**：设观察线法线归一化后逆时针旋转 90° 得单位切向量 `t`，则目标点 = 观察线给定点 + `target_offset · t`；**目标偏差** = 观察线交点相对目标点沿 `t` 的有符号距离（即求解目标偏差为 0 的发射角）。
- **求解流程**：在 `[-half_angle_deg, +half_angle_deg]` 内含端点等间距取 `sample_count` 个发射角逐一追迹求偏差；仅对**相邻有效采样**（追迹 `completed` 且与观察线有前向交点）偏差**异号**的区间做二分；二分过程中遇无效追迹或无前向交点即丢弃该区间，**不跨无效采样连根**。
- **停止条件**：|目标偏差| ≤ 1e-6，或角度区间宽 ≤ 1e-7 度；单个区间最多迭代 60 次。采样点 |偏差| ≤ 1e-6 时直接记为解。
- **返回**：按角度升序、按 1e-7 度去重后的解列表；无解时返回空列表 `[]`。

| 方法 | 路径 | 说明 |
| --- | --- | --- |
| POST | `/api/aim-solutions` | 扫描 + 异号区间二分求解命中目标偏移量的发射角，返回角度、出射方向、观察线交点与目标偏差 |

### 瞄准求解请求示例

```bash
curl -X POST http://localhost:8000/api/aim-solutions \
  -H "Content-Type: application/json" \
  -d '{
    "start": {"x": 0, "y": 0},
    "center_direction": {"x": 1, "y": 0},
    "half_angle_deg": 10,
    "sample_count": 21,
    "n_start": 1.0,
    "interfaces": [
      {"point": {"x": 2, "y": 0}, "normal": {"x": 1, "y": 0}, "n_next": 1.5},
      {"point": {"x": 4, "y": 0}, "normal": {"x": 1, "y": 0}, "n_next": 1.0}
    ],
    "observation_line": {"point": {"x": 10, "y": 0}, "normal": {"x": 1, "y": 0}},
    "target_offset": 1.0
  }'
```

### 瞄准求解响应说明

```json
[
  {
    "emission_angle_deg": 6.118255615234375,
    "exit_direction": {"x": 0.994304035740236, "y": 0.1065808824821775},
    "observation_point": {"x": 10.0, "y": 0.999999463241611},
    "target_deviation": -5.367583889581695e-07
  }
]
```

- 响应为解的**数组**（无解时为 `[]`），按 `emission_angle_deg` 升序排列。
- `emission_angle_deg`：相对中心方向的发射角（度，逆时针为正）。
- `exit_direction`：穿过所有界面后的出射段方向（归一化）；无界面时即发射方向。
- `observation_point`：出射段与观察线的前向交点。
- `target_deviation`：该解的目标偏差，|值| ≤ 1e-6（因角度区间宽 ≤ 1e-7 度而停止时可能略超）。

## 装调灵敏度（法线角度导数）

`POST /api/alignment-sensitivities` 接收**单条光线**（起点、方向、起始折射率）、按顺序排列的界面、**观察线**（线上一点 + 法线）与**正值扰动角** `perturbation_angle_deg`（度，> 0），定量评估每个界面的装调角度误差对观察线命中位置的影响。该接口为即时计算，结果不写入数据库。

- **切向量与基准偏移**：将观察线法线归一化后**逆时针旋转 90°** 得单位切向量 `t = (-ny, nx)`；先做一次基准追迹，取最终出射段与观察线的**前向**交点 `P`，基准偏移量 `offset = (P − 观察线给定点) · t`（沿 `t` 为正的有符号距离）。
- **逐界面法线扰动**：对第 `i` 个界面，**只将其法线**分别旋转 `+δ` 与 `−δ`（`δ = perturbation_angle_deg`，逆时针为正；界面点、折射率与其他界面参数保持不变）后重新追迹，各得到一次偏移量。
- **中心差分导数**：`derivative_per_deg = (offset(+δ) − offset(−δ)) / (2δ)`，单位为**长度/度**。
- **失效处理（不跨失效分支求差分）**：基准或任一扰动追迹发生**全反射**（`total_internal_reflection`）、**漏过界面**（`no_intersection`，平行或交点在反方向），或出射段**无法前向命中观察线**（`no_observation_intersection`，平行或交点在光线反方向）时，该次 `offset` 与 `observation_point` 为 `null`、`reason` 给出原因；只要 `+δ`/`−δ` 任一侧失效，该界面的 `derivative_per_deg` 即为 `null`。注意：基准失效不会跳过扰动计算，每个界面仍返回两次扰动的完整状态，只是导数恒为 `null`。
- 扰动结果中附 `perturbed_normal`（旋转后的单位法线）与 `perturbation_angle_deg`（`+δ` / `−δ`），便于核对。

| 方法 | 路径 | 说明 |
| --- | --- | --- |
| POST | `/api/alignment-sensitivities` | 基准偏移 + 逐界面 ±法线扰动追迹，返回两次扰动结果与中心差分导数（长度/度，不入库） |

### 装调灵敏度请求示例

```bash
curl -X POST http://localhost:8000/api/alignment-sensitivities \
  -H "Content-Type: application/json" \
  -d '{
    "start": {"x": 0, "y": 0},
    "direction": {"x": 1, "y": 0},
    "n_start": 1.0,
    "interfaces": [
      {"point": {"x": 2, "y": 0}, "normal": {"x": 1, "y": 0}, "n_next": 1.5},
      {"point": {"x": 4, "y": 0}, "normal": {"x": 1, "y": 0}, "n_next": 1.0}
    ],
    "observation_line": {"point": {"x": 10, "y": 0}, "normal": {"x": 1, "y": 0}},
    "perturbation_angle_deg": 0.1
  }'
```

### 装调灵敏度响应说明

```json
{
  "tangent": {"x": 0.0, "y": 1.0},
  "perturbation_angle_deg": 0.1,
  "baseline": {
    "status": "completed",
    "reason": null,
    "observation_point": {"x": 10.0, "y": 0.0},
    "offset": 0.0
  },
  "interfaces": [
    {
      "interface_index": 0,
      "plus": {
        "status": "completed",
        "reason": null,
        "observation_point": {"x": 10.0, "y": 0.00639954603029729},
        "offset": 0.00639954603029729,
        "perturbed_normal": {"x": 0.9999984769132877, "y": 0.0017453283658983088},
        "perturbation_angle_deg": 0.1
      },
      "minus": {
        "status": "completed",
        "reason": null,
        "observation_point": {"x": 10.0, "y": -0.00639954603029729},
        "offset": -0.00639954603029729,
        "perturbed_normal": {"x": 0.9999984769132877, "y": -0.0017453283658983088},
        "perturbation_angle_deg": -0.1
      },
      "derivative_per_deg": 0.06399546030297289
    }
  ]
}
```

- `tangent`：观察线单位切向量（法线逆时针旋转 90°）；`perturbation_angle_deg` 回显扰动角 δ。
- `baseline`：基准追迹状态 `status`（`completed` / `total_internal_reflection` / `no_intersection` / `no_observation_intersection`）、停止原因 `reason`、观察线前向交点 `observation_point` 与基准有符号偏移 `offset`；失效时后两者为 `null`。
- `interfaces[]` 与请求界面一一对应：`interface_index` 为序号，`plus` / `minus` 分别为法线旋转 `+δ` / `−δ` 后的追迹结果（字段同基准，另含扰动后的单位法线与实际扰动角）；`derivative_per_deg` 为中心差分导数（长度/度），任一侧或基准失效时为 `null`。

服务地址：API 根路径 http://localhost:8000 ，交互式文档（可在线试调该接口）http://localhost:8000/docs 。

## 项目结构

```
app/
  main.py       FastAPI 路由与生命周期（启动时建库）
  geometry.py   二维求交、斯涅尔折射、方向采样、光斑质心/RMS、菲涅耳 s/p 透过率与功率加权质心、瞄准求解（扫描+二分）、装调灵敏度（逐界面法线 ±扰动与中心差分）
  db.py         SQLite 持久化（traces / segments 两表）
  schemas.py    Pydantic 请求 / 响应模型
Dockerfile
docker-compose.yml
requirements.txt
```
