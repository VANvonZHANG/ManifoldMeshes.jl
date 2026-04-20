# sphere/ — 球面网格实现

LatLonGrid: 结构化经纬网格，全球覆盖 `[-90°, 90°] × [0°, 360°]`。

## 数据布局

- `nodes::Matrix{SVector{3,Float64}}` — `(nlat+1) × (nlon+1)`，nodes[1, :] = 南极
- `cell_volumes::Matrix{Float64}` — `nlat × nlon`，构造时预计算
- `cell_centroids::Matrix{SVector{3,Float64}}` — `nlat × nlon`，构造时预计算
- 线性索引: `cell_id = (ilat-1) * nlon + ilon`
- 节点索引: `node_id = (ilat-1) * (nlon+1) + ilon`
- 边编号: 水平 1..(nlat+1)×nlon，垂直 (nlat+1)×nlon+1..total

## 极点退化处理

球面网格的核心难点。北极/南极处多个节点退化为同一物理点。

**cell_volume 对角线退化:**
- 每个球面四边形沿 A-C 对角线 (SW→NE) 拆为 2 个球面三角形
- 若 A-C 退化 (A==C 或 B==D 反平行导致三角形跨半圆)，改用 B-D 对角线
- 若两条对角线都退化 (180° 极地单元，顶点退化为 2 个对跖点)，使用 lune 公式: `R² × (sin(lat₂) - sin(lat₁)) × |Δlon|`

**edge_outward_normal 极点坍缩:**
- 零长度边 (极点处重合端点) 返回 `zero(SVector{3,Float64})`

## 周期边界

经度方向通过取模实现东西环绕:
- `ilon_next = (ilon % nlon) + 1`
- `cell_cells` 东西向邻居自动环绕，南北向边界返回 0 (sentinel)

## l'Huilier 球面三角形面积

构造时热循环，用原始点积而非 `Manifolds.distance` 提升性能。数值稳定: clamp `acos` 参数到 `[-1, 1]`，`max(tan_half, 0.0)` 防止负数开方。

## 依赖 Manifolds.jl 的调用

- `Manifolds.mean(manifold, vertices)` — 计算 cell_centroid (Fréchet mean)
- `Manifolds.mid_point(manifold, n1, n2)` — edge_midpoint
- `Manifolds.project(manifold, base_point, vector)` — edge_outward_normal 投影到切空间
