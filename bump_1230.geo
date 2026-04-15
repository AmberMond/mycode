// 参数定义
L = 3.0;      // 通道长度
H = 1.0;      // 通道高度
BumpH = 0.1;  // 凸起高度 (10% 弦长)
Lc = 0.05;    // 网格特征长度 (越小网格越密)

// 1. 定义点 (Points)
Point(1) = {-1.5, 0, 0, Lc};  // 入口下角点
Point(2) = {-0.5, 0, 0, Lc};  // 凸起开始
Point(3) = {0, BumpH, 0, Lc}; // 凸起顶点
Point(4) = {0.5, 0, 0, Lc};   // 凸起结束
Point(5) = {1.5, 0, 0, Lc};   // 出口下角点
Point(6) = {1.5, H, 0, Lc};   // 出口上角点
Point(7) = {-1.5, H, 0, Lc};  // 入口上角点

// 2. 定义线 (Lines)
Line(1) = {1, 2};             // 下壁面左段
Spline(2) = {2, 3, 4};        // 凸起 (使用样条曲线)
Line(3) = {4, 5};             // 下壁面右段
Line(4) = {5, 6};             // 出口
Line(5) = {6, 7};             // 上壁面
Line(6) = {7, 1};             // 入口

// 3. 定义线环 (Curve Loop) 和 平面 (Plane Surface)
Curve Loop(1) = {1, 2, 3, 4, 5, 6};
Plane Surface(1) = {1};

// 4. 定义物理组 (Physical Groups) - **这对求解器至关重要**
// 求解器通过这些 ID 识别边界条件
Physical Curve("Inlet", 1) = {6};          // ID 1: 入口
Physical Curve("Outlet", 2) = {4};         // ID 2: 出口
Physical Curve("Wall", 3) = {1, 2, 3, 5};  // ID 3: 固壁 (上下)
Physical Surface("Domain", 100) = {1};     // ID 100: 内部流场