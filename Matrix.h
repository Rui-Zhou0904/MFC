#pragma once
#define M_PI 3.1415926
class CMatrix
{
public:
	CMatrix();
	CMatrix(int rowN, int colN);
	CMatrix(const CMatrix& m);
	~CMatrix(void);
private:
	double **dMatData;    //保存矩阵元素数据的二维数组
	int iRow;             //矩阵的行
	int iCol;             //矩阵的列
	// 基本操作
public:
	int Row() const { return iRow; };                  //返回行
	int Col() const { return iCol; };                  //返回列
	void SetSize(int rowN, int colN);                  //调整数组的大小，原有数据不变
	double& operator () (int row, int col);          //获取矩阵元素
	double operator () (int row, int col) const;     //重载获取矩阵元素函数，只有const对象能访问
	CMatrix& operator = (const CMatrix& m);
	friend CMatrix operator + (const CMatrix& m1, const CMatrix& m2);    //注意：友元函数并不是类自己的成员函数
	friend CMatrix operator - (const CMatrix& m1, const CMatrix& m2);
	friend CMatrix operator * (const CMatrix& m1, const CMatrix& m2);
	friend CMatrix operator * (const double& num, const CMatrix& m1);
	friend CMatrix operator * (const CMatrix& m1, const double& num);
	friend CMatrix operator ~ (const CMatrix& m);    //矩阵转置
	void Unit();      //生成单位矩阵
	// 矩阵操作
	void RowExtraction(CMatrix m, int row);       // 行提取
	void ColExtraction(CMatrix m, int col);       // 列提取
	int* RowOrder(int row);                       // 行由大到小排序
	int* ColOrder(int col);                       // 列由大到小排序
	int* MatrixMax();                             // 矩阵最大值
	// 二维空间运算
	void Vector2D(double X, double Y);                         // 生成二维行向量
	void NormVector2D();                                       // 标准化的二维行向量
	void NormVector3D();                                       // 标准化的三维行向量
	double Vector2DValue();                                    // 二维行向量长度
	double Vector3DValue();                                    // 三维行向量长度
	void ACRotate2D(double Phi);                               // 二维旋转矩阵
	void CRotate2D(double Phi);                                // 二维旋转矩阵
	double Clock2D(const CMatrix& m1);                         // 两个二维列向量的顺时针夹角，自身和输入都必须为二维列向量
	double AntiClock2D(const CMatrix& m1);                     // 两个二维列向量的逆时针夹角，自身和输入都必须为二维列向量
	void RegulateVector2(double Angle, double Radius);         // 给定二维向量，调整调整方向和长度
	// 航迹规划函数
	void Dubins(CMatrix Start, double SD, CMatrix End, double ED, double Rmin);
	void FormationDubins(CMatrix Start, double SD, CMatrix End, double ED, double Rmin, double UAVInitialFormationY, double FormationYLower, double FormationYUpper);
	// 沿航迹移动函数
	void PathMove(double Curvature, double Power);
	// 坐标转换
	void CarToGeo(double Lon0, double Lat0, double X, double Y);
	void GeoToCar(double Lon0, double Lat0, double Lon, double Lat);
};

