// pm_NASA_exact.geo
L_domain = 2.2;
H_domain = 1.0;
x_corner = 1.0;
angle_rad = Atan(0.321539 / 1.2); 

Point(1) = {0, 0, 0};
Point(2) = {x_corner, 0, 0};
Point(3) = {L_domain, -(L_domain - x_corner) * Tan(angle_rad), 0};
Point(4) = {L_domain, H_domain, 0};
Point(5) = {0, H_domain, 0};

// 定义边界线
Line(1) = {1, 2}; // 前缘无黏壁面 (从入口底端指向拐角)
Line(2) = {2, 3}; // 膨胀扇无黏壁面 (从拐角指向出口底端)
Line(3) = {3, 4}; // 出口界
Line(4) = {4, 5}; // 顶侧远场
Line(5) = {5, 1}; // 入口界

// -------------------------------------------------------------
// 重点：使用 Transfinite Curve 结合 Progression 进行几何级数布点
// -------------------------------------------------------------

// 1. 前缘壁面 (Line 1: 点1 -> 点2)
// 设置 15 个点，Progression < 1 代表间距越来越小（向拐角聚集）
Transfinite Curve {1} = 15 Using Progression 0.70;

// 2. 膨胀扇壁面 (Line 2: 点2 -> 点3)
// 设置 35 个点，Progression > 1 代表间距越来越大（渐渐远离拐角）
Transfinite Curve {2} = 35 Using Progression 1.15;

// 3. 入口与出口等粗大边界
Transfinite Curve {5} = 5;  // 入口只有5个点，和NASA图片吻合
Transfinite Curve {4} = 15; // 顶部远场，均匀分布
Transfinite Curve {3} = 15; // 出口，均匀分布

Curve Loop(1) = {1, 2, 3, 4, 5};
Plane Surface(1) = {1};

// 物理组
Physical Curve("Inlet", 1) = {5};
Physical Curve("Outlet", 2) = {3, 4};
Physical Curve("Inviscid_Wall", 3) = {1, 2};
Physical Surface("Domain", 4) = {1};

// -------------------------------------------------------------
// 算法设置
// -------------------------------------------------------------
// 取消背景场的干扰，依靠我们设定的边界点向内推进
Mesh.CharacteristicLengthExtendFromBoundary = 1;
// 使用标准的 Delaunay 算法 (Algorithm 5) 或 Frontal-Delaunay (Algorithm 6)
// 这两种算法在不加外部干涉时，会自动产生如图1所示的由密到疏的“辐射连线”
Mesh.Algorithm = 5; 