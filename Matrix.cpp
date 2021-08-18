#include "stdafx.h"
#include "Matrix.h"
#include "float.h"

CMatrix::CMatrix()
{
	iRow = 1; iCol = 1;
	dMatData = new double*[1];
	dMatData[0] = new double[1];
	dMatData[0][0] = 0;
}

CMatrix::CMatrix(int rowN, int colN)
{
	iRow = (rowN <= 0) ? 1 : rowN;
	iCol = (colN <= 0) ? 1 : colN;
	dMatData = new double*[iRow];
	for (int i = 0; i < iRow; i++)
	{
		dMatData[i] = new double[iCol];
		for (int j = 0; j < iCol; j++) { dMatData[i][j] = 0; }
	}
}

// 拷贝构造函数的作用：
//（1）以类对象作为函数参数传值调用时;
//（2）函数返回值为类对象；
//（3）用一个已定义的对象去初始化一个新对象时；
CMatrix::CMatrix(const CMatrix& m)
{
	iRow = m.Row(); iCol = m.Col();
	dMatData = new double*[iRow];
	for (int i = 0; i < iRow; i++)
	{
		dMatData[i] = new double[iCol];
		for (int j = 0; j < iCol; j++) { memcpy(dMatData[i], m.dMatData[i], sizeof(double)* iCol); }
	}
}

CMatrix::~CMatrix(void)
{
	for (int i = 0; i < iRow; i++) { delete[] dMatData[i]; }
	delete[] dMatData;
}

//调整矩阵大小，原有值不变
void CMatrix::SetSize(int rowN, int colN)
{
	if ((rowN <= 0) || (colN <= 0)) { return; }
	if ((rowN == iRow) && (colN == iCol)) { return; }
	double **rsData = new double*[rowN];
	for (int i = 0; i < rowN; i++)
	{
		rsData[i] = new double[colN];
		for (int j = 0; j < colN; j++)
		{
			rsData[i][j] = 0;
		}
	}
	int minRow = (iRow > rowN) ? rowN : iRow;
	int minCol = (iCol > colN) ? colN : iCol;
	int colSize = minCol * sizeof(double);
	for (int i = 0; i < minRow; i++)
	{
		memcpy(rsData[i], dMatData[i], colSize);
	}
	for (int i = 0; i < minRow; i++)
	{
		delete[] dMatData[i];
	}
	delete[] dMatData;
	dMatData = rsData;
	iRow = rowN; iCol = colN;
}

//返回数组元素（引用返回）
double& CMatrix::operator () (int row, int col)
{
	int i = (row < 0) ? 0 : row; i = (i >(iRow - 1)) ? (iRow - 1) : i;
	int j = (col < 0) ? 0 : col; j = (j >(iCol - 1)) ? (iCol - 1) : j;
	return dMatData[i][j];
}

//返回数组元素（重载）
double CMatrix::operator () (int row, int col) const
{
	int i = (row < 0) ? 0 : row; i = (i >(iRow - 1)) ? (iRow - 1) : i;
	int j = (col < 0) ? 0 : col; j = (j >(iCol - 1)) ? (iCol - 1) : j;
	return dMatData[i][j];
}

//重载赋值运算符=，当左右两边矩阵的大小不相等时，
//以右边的大小为基准,调整左边矩阵的大小
CMatrix &CMatrix::operator = (const CMatrix& m)
{
	if (iRow != m.Row() || iCol != m.Col())
	{
		SetSize(m.Row(), m.Col());
	}
	for (int i = 0; i < iRow; i++)
	{
		for (int j = 0; j < iCol; j++)
		{
			dMatData[i][j] = m(i, j);
		}
	}
	return *this;
}

//重载预算符+
CMatrix operator + (const CMatrix& m1, const CMatrix& m2)
{
	CMatrix CM1(m1); CMatrix CM2(m2);
	int rowN = max(CM1.Row(), CM2.Row());
	int colN = max(CM1.Col(), CM2.Col());
	CM1.SetSize(rowN, colN); CM2.SetSize(rowN, colN);
	CMatrix CM(rowN, colN);
	for (int i = 0; i < rowN; i++)
	{
		for (int j = 0; j < colN; j++)
		{
			CM.dMatData[i][j] = CM1(i, j) + CM2(i, j);
		}
	}
	return CM;
}

//重载预算符-
CMatrix operator - (const CMatrix& m1, const CMatrix& m2)
{
	CMatrix CM1(m1); CMatrix CM2(m2);
	int rowN = max(CM1.Row(), CM2.Row());
	int colN = max(CM1.Col(), CM2.Col());
	CM1.SetSize(rowN, colN); CM2.SetSize(rowN, colN);
	CMatrix CM(rowN, colN);
	for (int i = 0; i < rowN; i++)
	{
		for (int j = 0; j < colN; j++)
		{
			CM.dMatData[i][j] = CM1(i, j) - CM2(i, j);
		}
	}
	return CM;
}

//重载预算符*,两个矩阵相乘，m1的列要等于m2的行
CMatrix operator * (const CMatrix& m1, const CMatrix& m2)
{
	CMatrix CM1(m1); CMatrix CM2(m2);
	int CR = max(CM1.Col(), CM2.Row());
	CM1.SetSize(CM1.Row(), CR);
	CM2.SetSize(CR, CM2.Col());
	CMatrix CM(CM1.Row(), CM2.Col());
	for (int i = 0; i < CM1.Row(); i++)
	{
		for (int j = 0; j < CM2.Col(); j++)
		{
			for (int k = 0; k < CR; k++)
			{
				CM.dMatData[i][j] += CM1(i, k) * CM2(k, j);
			}
		}
	}
	return CM;
}

//重载预算符*,矩阵右乘一个数
CMatrix operator * (const CMatrix& m1, const double& num)
{
	CMatrix matTmp(m1.Row(), m1.Col());
	for (int i = 0; i < m1.Row(); i++)
	{
		for (int j = 0; j < m1.Col(); j++)
		{
			matTmp.dMatData[i][j] = m1(i, j) * num;
		}
	}
	return matTmp;
}

//重载预算符*,矩阵左乘一个数
CMatrix operator * (const double& num, const CMatrix& m1)
{
	CMatrix matTmp(m1.Row(), m1.Col());
	for (int i = 0; i < m1.Row(); i++)
	{
		for (int j = 0; j < m1.Col(); j++)
		{
			matTmp.dMatData[i][j] = m1(i, j) * num;
		}
	}
	return matTmp;
}

//矩阵转置
CMatrix operator ~ (const CMatrix& m)
{
	CMatrix matTmp(m.Col(), m.Row());
	for (int i = 0; i < m.Row(); i++)
	{
		for (int j = 0; j < m.Col(); j++)
		{
			matTmp.dMatData[i][j] = m(i, j);
		}
	}
	return matTmp;
}

//单位化矩阵
void CMatrix::Unit()
{
	for (int i = 0; i < iRow; i++)
	{
		for (int j = 0; j < iCol; j++)
		{
			dMatData[i][j] = (i == j) ? 1 : 0;
		}
	}
}

// 行提取
void CMatrix::RowExtraction(CMatrix m, int row)
{
	SetSize(1, m.Col());
	int i = (row < 0) ? 0 : row; i = (i >(m.Row() - 1)) ? (m.Row() - 1) : i;
	for (int j = 0; j < m.Col(); j++)
	{
		dMatData[0][j] = m(i, j);
	}
}

// 列提取
void CMatrix::ColExtraction(CMatrix m, int col)
{
	SetSize(m.Row(), 1);
	int j = (col < 0) ? 0 : col; j = (j >(m.Col() - 1)) ? (m.Col() - 1) : j;
	for (int i = 0; i < m.Row(); i++)
	{
		dMatData[i][0] = m(i, j);
	}
}

// 行前N个最大
int* CMatrix::RowOrder(int RowNum)
{
	int R = (RowNum < 0) ? 0 : RowNum; R = (R > iRow) ? iRow : R;
	static int* Order = new int[iCol];
	double* TheRow = new double[iCol];
	for (int i = 0; i < iCol; i++)
	{
		TheRow[i] = dMatData[R][i];
		Order[i] = i;
	}
	for (int i = 0; i < iCol - 1; i++)
	{
		for (int j = 0; j < iCol - i - 1; j++)
		{
			if (TheRow[j] < TheRow[j + 1])
			{
				double tempD = TheRow[j + 1]; TheRow[j + 1] = TheRow[j]; TheRow[j] = tempD;
				int tempI = Order[j + 1]; Order[j + 1] = Order[j]; Order[j] = tempI;
			}
		}
	}
	return Order;
}

// 列前N个最大
int* CMatrix::ColOrder(int ColNum)
{
	int C = (ColNum < 0) ? 0 : ColNum; C = (C > iCol) ? iCol : C;
	static int* Order = new int[iRow];
	double* TheCol = new double[iRow];
	for (int i = 0; i < iRow; i++)
	{
		TheCol[i] = dMatData[i][C];
		Order[i] = i;
	}
	for (int i = 0; i < iRow - 1; i++)
	{
		for (int j = 0; j < iRow - i - 1; j++)
		{
			if (TheCol[j] < TheCol[j + 1])
			{
				double tempD = TheCol[j + 1]; TheCol[j + 1] = TheCol[j]; TheCol[j] = tempD;
				int tempI = Order[j + 1]; Order[j + 1] = Order[j]; Order[j] = tempI;
			}
		}
	}
	return Order;
}

int* CMatrix::MatrixMax()
{
	static int TheMax[2];
	double tempMax = -1 * DBL_MAX;
	for (int i = 0; i < iRow; i++)
	{
		for (int j = 0; j < iCol; j++)
		{
			if (dMatData[i][j] > tempMax)
			{
				tempMax = dMatData[i][j];
				TheMax[0] = i; TheMax[1] = j;
			}
		}
	}
	return TheMax;
}

// 生成二维行向量
void CMatrix::Vector2D(double X, double Y)
{
	SetSize(1, 2);
	dMatData[0][0] = X; dMatData[0][1] = Y;
}

// 标准化的二维行向量
void CMatrix::NormVector2D()
{
	SetSize(1, 2);
	double norm = sqrt(dMatData[0][0] * dMatData[0][0] + dMatData[0][1] * dMatData[0][1]);
	dMatData[0][0] = dMatData[0][0] / norm; dMatData[0][1] = dMatData[0][1] / norm;
}

// 标准化的三维行向量
void CMatrix::NormVector3D()
{
	SetSize(1, 3);
	double norm = sqrt(dMatData[0][0] * dMatData[0][0] + dMatData[0][1] * dMatData[0][1] + dMatData[0][2] * dMatData[0][2]);
	dMatData[0][0] = dMatData[0][0] / norm; dMatData[0][1] = dMatData[0][1] / norm; dMatData[0][2] = dMatData[0][2] / norm;
}

double CMatrix::Vector2DValue()
{
	SetSize(1, 2);
	return sqrt(dMatData[0][0] * dMatData[0][0] + dMatData[0][1] * dMatData[0][1]);
}

double CMatrix::Vector3DValue()
{
	SetSize(1, 3);
	return sqrt(dMatData[0][0] * dMatData[0][0] + dMatData[0][1] * dMatData[0][1] + dMatData[0][2] * dMatData[0][2]);
}

// 二维旋转矩阵
void CMatrix::ACRotate2D(double Phi)
{
	SetSize(2, 2);
	dMatData[0][0] = cos(Phi);  dMatData[0][1] = sin(Phi);
	dMatData[1][0] = -sin(Phi); dMatData[1][1] = cos(Phi);
}

// 二维旋转矩阵
void CMatrix::CRotate2D(double Phi)
{
	SetSize(2, 2);
	dMatData[0][0] = cos(Phi);  dMatData[0][1] = -sin(Phi);
	dMatData[1][0] = sin(Phi); dMatData[1][1] = cos(Phi);
}

// 两个二维列向量的顺时针夹角
double CMatrix::Clock2D(const CMatrix& m1)
{
	SetSize(1, 2); CMatrix CM1(m1); CM1.SetSize(1, 2);
	double Cal = (dMatData[0][0] * CM1(0, 0) + dMatData[0][1] * CM1(0, 1)) / \
		(sqrt(dMatData[0][0] * dMatData[0][0] + dMatData[0][1] * dMatData[0][1]) * sqrt(CM1(0, 0) * CM1(0, 0) + CM1(0, 1) * CM1(0, 1)));
	Cal = (Cal > 1) ? 1 : Cal; Cal = (Cal < -1) ? (-1) : Cal;
	double Angle = ((dMatData[0][0] * CM1(0, 1) - dMatData[0][1] * CM1(0, 0)) >= 0) ? (-acos(Cal)) : acos(Cal);
	if (Angle < 0) { Angle += (2 * M_PI); }
	return Angle;
}

// 两个二维列向量的逆时针夹角
double CMatrix::AntiClock2D(const CMatrix& m1)
{
	SetSize(1, 2); CMatrix CM1(m1); CM1.SetSize(1, 2);
	double Cal = (dMatData[0][0] * CM1(0, 0) + dMatData[0][1] * CM1(0, 1)) / \
		(sqrt(dMatData[0][0] * dMatData[0][0] + dMatData[0][1] * dMatData[0][1]) * sqrt(CM1(0, 0) * CM1(0, 0) + CM1(0, 1) * CM1(0, 1)));
	Cal = (Cal > 1) ? 1 : Cal; Cal = (Cal < -1) ? (-1) : Cal;
	double Angle = ((dMatData[0][0] * CM1(0, 1) - dMatData[0][1] * CM1(0, 0)) >= 0) ? acos(Cal) : (-acos(Cal));
	if (Angle < 0) { Angle += (2 * M_PI); }
	return Angle;
}

// 给定二维向量，调整调整方向和长度
void CMatrix::RegulateVector2(double Angle, double Radius)
{
	double x = dMatData[0][0];
	double y = dMatData[0][1];
	double TheSin = sin(M_PI * Angle / 180);
	double TheCos = cos(M_PI * Angle / 180);
	double newX = x * TheCos + -y * TheSin;
	double newY = x * TheSin + y * TheCos;
	CMatrix TempVector(1, 2); TempVector(0, 0) = newX; TempVector(0, 1) = newY;
	double K = Radius / TempVector.Vector2DValue();
	dMatData[0][0] = newX * K; dMatData[0][1] = newY * K;
}

// 二维Dubins
void CMatrix::Dubins(CMatrix Start, double SD, CMatrix End, double ED, double Rmin)
{
	SetSize(3, 2);
	CMatrix RotationMatrixA(2, 2), RotationMatrixB(2, 2), Arrow(1, 2), tempVector(1, 2);
	double RminD = (double)Rmin; int DubinsType; double DubinsLength;
	double AAX, AAY, AAMR;
	//double AAM;
	double LSL_L; double LSLCurvatureVector[3], LSLLengthVector[3];
	double LSR_L; double LSRCurvatureVector[3], LSRLengthVector[3]; bool IfLSR;
	double RSL_L; double RSLCurvatureVector[3], RSLLengthVector[3]; bool IfRSL;
	double RSR_L; double RSRCurvatureVector[3], RSRLengthVector[3];
	//double LRL_L; double LRLCurvatureVector[3], LRLLengthVector[3]; bool IfLRL;
	//double RLR_L; double RLRCurvatureVector[3], RLRLengthVector[3]; bool IfRLR;
	// 建立起点、终点
	CMatrix StartXY(Start); StartXY.SetSize(1, 2);
	CMatrix EndXY(End); EndXY.SetSize(1, 2);
	// 建立四个adjacent circle
	CMatrix AS(1, 2); AS.Vector2D(cos(M_PI * SD / 180.0), sin(M_PI * SD / 180.0)); AS.NormVector2D();
	CMatrix AE(1, 2); AE.Vector2D(cos(M_PI * ED / 180.0), sin(M_PI * ED / 180.0)); AE.NormVector2D();
	RotationMatrixA.ACRotate2D(M_PI / 2);
	CMatrix CXL(1, 2); CXL = StartXY + RminD * AS * RotationMatrixA;
	CMatrix CXR(1, 2); CXR = StartXY - RminD * AS * RotationMatrixA;
	CMatrix CYL(1, 2); CYL = EndXY + RminD * AE * RotationMatrixA;
	CMatrix CYR(1, 2); CYR = EndXY - RminD * AE * RotationMatrixA;
	CMatrix DirectionXL(1, 2); DirectionXL = StartXY - CXL;
	CMatrix DirectionXR(1, 2); DirectionXR = StartXY - CXR;
	CMatrix DirectionYL(1, 2); DirectionYL = EndXY - CYL;
	CMatrix DirectionYR(1, 2); DirectionYR = EndXY - CYR;
	// 求LSL路线
	RotationMatrixA.ACRotate2D(-M_PI / 2);
	tempVector = CYL - CXL;
	AAMR = tempVector.Vector2DValue();
	tempVector.NormVector2D();
	Arrow = RminD * tempVector * RotationMatrixA;
	AAX = DirectionXL.AntiClock2D(Arrow);
	AAY = Arrow.AntiClock2D(DirectionYL);
	if (AAX > (2 * M_PI - 0.05)) { AAX = 0.0; }    // 防止大圆圈出现
	if (AAY > (2 * M_PI - 0.05)) { AAY = 0.0; }    // 防止大圆圈出现
	LSLLengthVector[0] = AAX * RminD; LSLLengthVector[1] = AAMR; LSLLengthVector[2] = AAY * RminD;
	LSLCurvatureVector[0] = 1 / RminD; LSLCurvatureVector[1] = 0; LSLCurvatureVector[2] = 1 / RminD;
	LSL_L = LSLLengthVector[0] + LSLLengthVector[1] + LSLLengthVector[2];
	// 求LSR路线
	tempVector = CYR - CXL;
	IfLSR = (tempVector.Vector2DValue() >= (2 * RminD));
	if (IfLSR)
	{
		double Angle = acos(2 * RminD / tempVector.Vector2DValue());
		tempVector.NormVector2D();
		RotationMatrixA.ACRotate2D(-Angle);
		CMatrix ArrowX(1, 2); ArrowX = RminD * tempVector * RotationMatrixA;
		CMatrix ArrowY(1, 2); ArrowY = ArrowY - ArrowX;
		CMatrix TangentYR_XL(1, 2); TangentYR_XL = CYR - CXL + ArrowY - ArrowX;
		AAMR = TangentYR_XL.Vector2DValue();
		AAX = DirectionXL.AntiClock2D(ArrowX);
		AAY = 2 * M_PI - ArrowY.AntiClock2D(DirectionYR);
		if (AAX > (2 * M_PI - 0.05)) { AAX = 0.0; }    // 防止大圆圈出现
		if (AAY > (2 * M_PI - 0.05)) { AAY = 0.0; }    // 防止大圆圈出现
		LSRLengthVector[0] = AAX * RminD; LSRLengthVector[1] = AAMR; LSRLengthVector[2] = AAY * RminD;
		LSRCurvatureVector[0] = 1 / RminD; LSRCurvatureVector[1] = 0; LSRCurvatureVector[2] = -1 / RminD;
		LSR_L = LSRLengthVector[0] + LSRLengthVector[1] + LSRLengthVector[2];
	}
	// 求RSL路线
	tempVector = CYL - CXR;
	IfRSL = (tempVector.Vector2DValue() >= (2 * RminD));
	if (IfRSL)
	{
		double Angle = acos(2 * RminD / tempVector.Vector2DValue());
		tempVector.NormVector2D();
		RotationMatrixA.ACRotate2D(Angle);
		CMatrix ArrowX(1, 2); ArrowX = RminD * tempVector * RotationMatrixA;
		CMatrix ArrowY(1, 2); ArrowY = ArrowY - ArrowX;
		CMatrix TangentYL_XR(1, 2); TangentYL_XR = CYL - CXR + ArrowY - ArrowX;
		AAMR = TangentYL_XR.Vector2DValue();
		AAX = 2 * M_PI - DirectionXR.AntiClock2D(ArrowX);
		AAY = ArrowY.AntiClock2D(DirectionYL);
		if (AAX > (2 * M_PI - 0.05)) { AAX = 0.0; }    // 防止大圆圈出现
		if (AAY > (2 * M_PI - 0.05)) { AAY = 0.0; }    // 防止大圆圈出现
		RSLLengthVector[0] = AAX * RminD; RSLLengthVector[1] = AAMR; RSLLengthVector[2] = AAY * RminD;
		RSLCurvatureVector[0] = -1 / RminD; RSLCurvatureVector[1] = 0; RSLCurvatureVector[2] = 1 / RminD;
		RSL_L = RSLLengthVector[0] + RSLLengthVector[1] + RSLLengthVector[2];
	}
	// 求RSR路线
	RotationMatrixA.ACRotate2D(M_PI / 2);
	tempVector = CYR - CXR;
	AAMR = tempVector.Vector2DValue();
	tempVector.NormVector2D();
	Arrow = RminD * tempVector * RotationMatrixA;
	AAX = 2 * M_PI - DirectionXR.AntiClock2D(Arrow);
	AAY = 2 * M_PI - Arrow.AntiClock2D(DirectionYR);
	if (AAX > (2 * M_PI - 0.05)) { AAX = 0.0; }    // 防止大圆圈出现
	if (AAY > (2 * M_PI - 0.05)) { AAY = 0.0; }    // 防止大圆圈出现
	RSRLengthVector[0] = AAX * RminD; RSRLengthVector[1] = AAMR; RSRLengthVector[2] = AAY * RminD;
	RSRCurvatureVector[0] = -1 / RminD; RSRCurvatureVector[1] = 0; RSRCurvatureVector[2] = -1 / RminD;
	RSR_L = RSRLengthVector[0] + RSRLengthVector[1] + RSRLengthVector[2];
	// 求LRL路线
	//tempVector = CYL - CXL;
	//IfLRL = (tempVector.Vector2DValue() <= (4 * RminD));
	//if(IfLRL)
	//{
	//	double Armx = sqrt((2 * RminD) * (2 * RminD) - (tempVector.Vector2DValue() / 2) * (tempVector.Vector2DValue() / 2));
	//	double Angle = asin(Armx / (2 * RminD));
	//	RotationMatrixA.ACRotate2D(Angle); RotationMatrixB.ACRotate2D(-Angle);
	//	tempVector.NormVector2D();
	//	CMatrix ArrowX1(1, 2); ArrowX1 = RminD * tempVector * RotationMatrixA;
	//	CMatrix ArrowY1(1, 2); ArrowY1 = ArrowY1 - RminD * tempVector * RotationMatrixB;
	//	CMatrix ArrowX2(1, 2); ArrowX2 = RminD * tempVector * RotationMatrixB;
	//	CMatrix ArrowY2(1, 2); ArrowY2 = ArrowY2 - RminD * tempVector * RotationMatrixA;
	//	double AAX1 = DirectionXL.AntiClock2D(ArrowX1);
	//	double AAM1 = 2 * M_PI - ArrowX1.AntiClock2D(ArrowY1);
	//	double AAY1 = ArrowY1.AntiClock2D(DirectionYL);
	//	double AAX2 = DirectionXL.AntiClock2D(ArrowX2);
	//	double AAM2 = 2 * M_PI - ArrowX2.AntiClock2D(ArrowY2);
	//	double AAY2 = ArrowY2.AntiClock2D(DirectionYL);
	//	if((AAX1 + AAM1 + AAY1) <= (AAX2 + AAM2 + AAY2)) {AAX = AAX1; AAY = AAY1; AAM = AAM1;}
	//	else {AAX = AAX2; AAY = AAY2; AAM = AAM2;}
	//	LRLLengthVector[0] = AAX * RminD; LRLLengthVector[1] = AAM * RminD; LRLLengthVector[2] = AAY * RminD;
	//	LRLCurvatureVector[0] = 1 / RminD; LRLCurvatureVector[1] = -1 / RminD; LRLCurvatureVector[2] = 1 / RminD;
	//	LRL_L = LRLLengthVector[0] + LRLLengthVector[1] + LRLLengthVector[2];
	//}
	// 求RLR路线
	//tempVector = CYR - CXR;
	//IfRLR = (tempVector.Vector2DValue() <= (4 * RminD));
	//if(IfRLR)
	//{
	//	double Armx = sqrt((2 * RminD) * (2 * RminD) - (tempVector.Vector2DValue() / 2) * (tempVector.Vector2DValue() / 2));
	//	double Angle = asin(Armx / (2 * RminD));
	//	RotationMatrixA.ACRotate2D(Angle); RotationMatrixB.ACRotate2D(-Angle);
	//	tempVector.NormVector2D();
	//	CMatrix ArrowX1(1, 2); ArrowX1 = RminD * tempVector * RotationMatrixA;
	//	CMatrix ArrowY1(1, 2); ArrowY1 = ArrowY1 - RminD * tempVector * RotationMatrixB;
	//	CMatrix ArrowX2(1, 2); ArrowX2 = RminD * tempVector * RotationMatrixB;
	//	CMatrix ArrowY2(1, 2); ArrowY2 = ArrowY2 - RminD * tempVector * RotationMatrixA;
	//	double AAX1 = 2 * M_PI - DirectionXR.AntiClock2D(ArrowX1);
	//	double AAM1 = ArrowX1.AntiClock2D(ArrowY1);
	//	double AAY1 = 2 * M_PI - ArrowY1.AntiClock2D(DirectionYR);
	//	double AAX2 = 2 * M_PI - DirectionXR.AntiClock2D(ArrowX2);
	//	double AAM2 = ArrowX2.AntiClock2D(ArrowY2);
	//	double AAY2 = 2 * M_PI - ArrowY2.AntiClock2D(DirectionYR);
	//	if((AAX1 + AAM1 + AAY1) <= (AAX2 + AAM2 + AAY2)) {AAX = AAX1; AAY = AAY1; AAM = AAM1;}
	//	else {AAX = AAX2; AAY = AAY2; AAM = AAM2;}
	//	RLRLengthVector[0] = AAX * RminD; RLRLengthVector[1] = AAM * RminD; RLRLengthVector[2] = AAY * RminD;
	//	RLRCurvatureVector[0] = -1 / RminD; RLRCurvatureVector[1] = 1 / RminD; RLRCurvatureVector[2] = -1 / RminD;
	//	RLR_L = RLRLengthVector[0] + RLRLengthVector[1] + RLRLengthVector[2];
	//}
	// 求Dubins路线
	DubinsType = (LSL_L <= RSR_L) ? 1 : 4; DubinsLength = min(LSL_L, RSR_L);
	if (IfLSR) { if (LSR_L < DubinsLength) { DubinsLength = LSR_L; DubinsType = 2; } }
	if (IfRSL) { if (RSL_L < DubinsLength) { DubinsLength = RSL_L; DubinsType = 3; } }
	//if(IfLRL) {if(LRL_L < DubinsLength) {DubinsLength = LRL_L; DubinsType = 5;}}
	//if(IfRLR) {if(RLR_L < DubinsLength) {DubinsLength = RLR_L; DubinsType = 6;}}
	switch (DubinsType)
	{
	case 1:
		dMatData[0][0] = LSLLengthVector[0];    dMatData[1][0] = LSLLengthVector[1];    dMatData[2][0] = LSLLengthVector[2];
		dMatData[0][1] = LSLCurvatureVector[0]; dMatData[1][1] = LSLCurvatureVector[1]; dMatData[2][1] = LSLCurvatureVector[2];
		break;
	case 2:
		dMatData[0][0] = LSRLengthVector[0];    dMatData[1][0] = LSRLengthVector[1];    dMatData[2][0] = LSRLengthVector[2];
		dMatData[0][1] = LSRCurvatureVector[0]; dMatData[1][1] = LSRCurvatureVector[1]; dMatData[2][1] = LSRCurvatureVector[2];
		break;
	case 3:
		dMatData[0][0] = RSLLengthVector[0];    dMatData[1][0] = RSLLengthVector[1];    dMatData[2][0] = RSLLengthVector[2];
		dMatData[0][1] = RSLCurvatureVector[0]; dMatData[1][1] = RSLCurvatureVector[1]; dMatData[2][1] = RSLCurvatureVector[2];
		break;
	case 4:
		dMatData[0][0] = RSRLengthVector[0];    dMatData[1][0] = RSRLengthVector[1];    dMatData[2][0] = RSRLengthVector[2];
		dMatData[0][1] = RSRCurvatureVector[0]; dMatData[1][1] = RSRCurvatureVector[1]; dMatData[2][1] = RSRCurvatureVector[2];
		break;
		//case 5:
		//	dMatData[0][0] = LRLLengthVector[0];    dMatData[1][0] = LRLLengthVector[1];    dMatData[2][0] = LRLLengthVector[2];
		//	dMatData[0][1] = LRLCurvatureVector[0]; dMatData[1][1] = LRLCurvatureVector[1]; dMatData[2][1] = LRLCurvatureVector[2];
		//	break;
		//case 6:
		//	dMatData[0][0] = RLRLengthVector[0];    dMatData[1][0] = RLRLengthVector[1];    dMatData[2][0] = RLRLengthVector[2];
		//	dMatData[0][1] = RLRCurvatureVector[0]; dMatData[1][1] = RLRCurvatureVector[1]; dMatData[2][1] = RLRCurvatureVector[2];
		//	break;
	}
}

// 编队二维Dubins
void CMatrix::FormationDubins(CMatrix Start, double SD, CMatrix End, double ED, double Rmin, double UAVInitialFormationY, double FormationYLower, double FormationYUpper)
{
	SetSize(3, 2);
	CMatrix RotationMatrixA(2, 2), RotationMatrixB(2, 2), Arrow(1, 2), tempVector(1, 2);
	double RminDL = (double)Rmin + (UAVInitialFormationY - FormationYLower); double RminDR = (double)Rmin + (FormationYUpper - UAVInitialFormationY);
	int DubinsType; double DubinsLength;
	double AAX, AAY, AAMR;
	double LSL_L; double LSLCurvatureVector[3], LSLLengthVector[3];
	double LSR_L; double LSRCurvatureVector[3], LSRLengthVector[3]; bool IfLSR;
	double RSL_L; double RSLCurvatureVector[3], RSLLengthVector[3]; bool IfRSL;
	double RSR_L; double RSRCurvatureVector[3], RSRLengthVector[3];
	// 建立起点、终点
	CMatrix StartXY(Start); StartXY.SetSize(1, 2);
	CMatrix EndXY(End); EndXY.SetSize(1, 2);
	// 建立四个adjacent circle
	CMatrix AS(1, 2); AS.Vector2D(cos(M_PI * SD / 180.0), sin(M_PI * SD / 180.0)); AS.NormVector2D();
	CMatrix AE(1, 2); AE.Vector2D(cos(M_PI * ED / 180.0), sin(M_PI * ED / 180.0)); AE.NormVector2D();
	RotationMatrixA.ACRotate2D(M_PI / 2);
	CMatrix CXL(1, 2); CXL = StartXY + RminDL * AS * RotationMatrixA;
	CMatrix CXR(1, 2); CXR = StartXY - RminDR * AS * RotationMatrixA;
	CMatrix CYL(1, 2); CYL = EndXY + RminDL * AE * RotationMatrixA;
	CMatrix CYR(1, 2); CYR = EndXY - RminDR * AE * RotationMatrixA;
	CMatrix DirectionXL(1, 2); DirectionXL = StartXY - CXL;
	CMatrix DirectionXR(1, 2); DirectionXR = StartXY - CXR;
	CMatrix DirectionYL(1, 2); DirectionYL = EndXY - CYL;
	CMatrix DirectionYR(1, 2); DirectionYR = EndXY - CYR;
	// 求LSL路线
	RotationMatrixA.ACRotate2D(-M_PI / 2);
	tempVector = CYL - CXL;
	AAMR = tempVector.Vector2DValue();
	tempVector.NormVector2D();
	Arrow = RminDL * tempVector * RotationMatrixA;
	AAX = DirectionXL.AntiClock2D(Arrow);
	AAY = Arrow.AntiClock2D(DirectionYL);
	if (AAX > (2 * M_PI - 0.05)) { AAX = 0.0; }    // 防止大圆圈出现
	if (AAY > (2 * M_PI - 0.05)) { AAY = 0.0; }    // 防止大圆圈出现
	LSLLengthVector[0] = AAX * RminDL; LSLLengthVector[1] = AAMR; LSLLengthVector[2] = AAY * RminDL;
	LSLCurvatureVector[0] = 1 / RminDL; LSLCurvatureVector[1] = 0; LSLCurvatureVector[2] = 1 / RminDL;
	LSL_L = LSLLengthVector[0] + LSLLengthVector[1] + LSLLengthVector[2];
	// 求LSR路线
	tempVector = CYR - CXL;
	IfLSR = (tempVector.Vector2DValue() >= (RminDL + RminDR));
	if (IfLSR)
	{
		double Angle = acos((RminDL + RminDR) / tempVector.Vector2DValue());
		tempVector.NormVector2D();
		RotationMatrixA.ACRotate2D(-Angle);
		CMatrix ArrowX(1, 2); ArrowX = RminDL * tempVector * RotationMatrixA;
		CMatrix ArrowY(1, 2); ArrowY = (RminDR / RminDL) * (ArrowY - ArrowX);
		CMatrix TangentYR_XL(1, 2); TangentYR_XL = CYR - CXL + ArrowY - ArrowX;
		AAMR = TangentYR_XL.Vector2DValue();
		AAX = DirectionXL.AntiClock2D(ArrowX);
		AAY = 2 * M_PI - ArrowY.AntiClock2D(DirectionYR);
		if (AAX > (2 * M_PI - 0.05)) { AAX = 0.0; }    // 防止大圆圈出现
		if (AAY > (2 * M_PI - 0.05)) { AAY = 0.0; }    // 防止大圆圈出现
		LSRLengthVector[0] = AAX * RminDL; LSRLengthVector[1] = AAMR; LSRLengthVector[2] = AAY * RminDR;
		LSRCurvatureVector[0] = 1 / RminDL; LSRCurvatureVector[1] = 0; LSRCurvatureVector[2] = -1 / RminDR;
		LSR_L = LSRLengthVector[0] + LSRLengthVector[1] + LSRLengthVector[2];
	}
	// 求RSL路线
	tempVector = CYL - CXR;
	IfRSL = (tempVector.Vector2DValue() >= (RminDL + RminDR));
	if (IfRSL)
	{
		double Angle = acos((RminDL + RminDR) / tempVector.Vector2DValue());
		tempVector.NormVector2D();
		RotationMatrixA.ACRotate2D(Angle);
		CMatrix ArrowX(1, 2); ArrowX = RminDR * tempVector * RotationMatrixA;
		CMatrix ArrowY(1, 2); ArrowY = (RminDL / RminDR) * (ArrowY - ArrowX);
		CMatrix TangentYL_XR(1, 2); TangentYL_XR = CYL - CXR + ArrowY - ArrowX;
		AAMR = TangentYL_XR.Vector2DValue();
		AAX = 2 * M_PI - DirectionXR.AntiClock2D(ArrowX);
		AAY = ArrowY.AntiClock2D(DirectionYL);
		if (AAX > (2 * M_PI - 0.05)) { AAX = 0.0; }    // 防止大圆圈出现
		if (AAY > (2 * M_PI - 0.05)) { AAY = 0.0; }    // 防止大圆圈出现
		RSLLengthVector[0] = AAX * RminDR; RSLLengthVector[1] = AAMR; RSLLengthVector[2] = AAY * RminDL;
		RSLCurvatureVector[0] = -1 / RminDR; RSLCurvatureVector[1] = 0; RSLCurvatureVector[2] = 1 / RminDL;
		RSL_L = RSLLengthVector[0] + RSLLengthVector[1] + RSLLengthVector[2];
	}
	// 求RSR路线
	RotationMatrixA.ACRotate2D(M_PI / 2);
	tempVector = CYR - CXR;
	AAMR = tempVector.Vector2DValue();
	tempVector.NormVector2D();
	Arrow = RminDR * tempVector * RotationMatrixA;
	AAX = 2 * M_PI - DirectionXR.AntiClock2D(Arrow);
	AAY = 2 * M_PI - Arrow.AntiClock2D(DirectionYR);
	if (AAX > (2 * M_PI - 0.05)) { AAX = 0.0; }    // 防止大圆圈出现
	if (AAY > (2 * M_PI - 0.05)) { AAY = 0.0; }    // 防止大圆圈出现
	RSRLengthVector[0] = AAX * RminDR; RSRLengthVector[1] = AAMR; RSRLengthVector[2] = AAY * RminDR;
	RSRCurvatureVector[0] = -1 / RminDR; RSRCurvatureVector[1] = 0; RSRCurvatureVector[2] = -1 / RminDR;
	RSR_L = RSRLengthVector[0] + RSRLengthVector[1] + RSRLengthVector[2];
	// 求Dubins路线
	DubinsType = (LSL_L <= RSR_L) ? 1 : 4; DubinsLength = min(LSL_L, RSR_L);
	if (IfLSR) { if (LSR_L < DubinsLength) { DubinsLength = LSR_L; DubinsType = 2; } }
	if (IfRSL) { if (RSL_L < DubinsLength) { DubinsLength = RSL_L; DubinsType = 3; } }
	switch (DubinsType)
	{
	case 1:
		dMatData[0][0] = LSLLengthVector[0];    dMatData[1][0] = LSLLengthVector[1];    dMatData[2][0] = LSLLengthVector[2];
		dMatData[0][1] = LSLCurvatureVector[0]; dMatData[1][1] = LSLCurvatureVector[1]; dMatData[2][1] = LSLCurvatureVector[2];
		break;
	case 2:
		dMatData[0][0] = LSRLengthVector[0];    dMatData[1][0] = LSRLengthVector[1];    dMatData[2][0] = LSRLengthVector[2];
		dMatData[0][1] = LSRCurvatureVector[0]; dMatData[1][1] = LSRCurvatureVector[1]; dMatData[2][1] = LSRCurvatureVector[2];
		break;
	case 3:
		dMatData[0][0] = RSLLengthVector[0];    dMatData[1][0] = RSLLengthVector[1];    dMatData[2][0] = RSLLengthVector[2];
		dMatData[0][1] = RSLCurvatureVector[0]; dMatData[1][1] = RSLCurvatureVector[1]; dMatData[2][1] = RSLCurvatureVector[2];
		break;
	case 4:
		dMatData[0][0] = RSRLengthVector[0];    dMatData[1][0] = RSRLengthVector[1];    dMatData[2][0] = RSRLengthVector[2];
		dMatData[0][1] = RSRCurvatureVector[0]; dMatData[1][1] = RSRCurvatureVector[1]; dMatData[2][1] = RSRCurvatureVector[2];
		break;
	}
}

// 沿航迹移动函数
void CMatrix::PathMove(double Curvature, double Power)
{
	SetSize(1, 3);
	CMatrix TempPoint(1, 2); TempPoint(0, 0) = dMatData[0][0]; TempPoint(0, 1) = dMatData[0][1];
	double TempDirection = dMatData[0][2];
	CMatrix TempVector(1, 2);
	CMatrix AS(1, 2); AS.Vector2D(cos(M_PI * TempDirection / 180.0), sin(M_PI * TempDirection / 180.0)); AS.NormVector2D();
	if (Curvature == 0) { TempVector = Power * AS; }
	else
	{
		CMatrix RotationMatrixA(2, 2); RotationMatrixA.ACRotate2D(M_PI / 2);
		CMatrix CX(1, 2); CX = (1 / Curvature) * AS * RotationMatrixA;
		CMatrix DirectionX(1, 2); DirectionX = (-1 / Curvature) * AS * RotationMatrixA;
		RotationMatrixA.ACRotate2D(Power * Curvature);
		TempVector = CX + DirectionX * RotationMatrixA;
		TempDirection += (180.0 * Power * Curvature / M_PI);
	}
	dMatData[0][0] += TempVector(0, 0); dMatData[0][1] += TempVector(0, 1); dMatData[0][2] = TempDirection;
}

void CMatrix::CarToGeo(double Lon0, double Lat0, double X, double Y)
{
	SetSize(1, 2);
	double tempDouble1 = 180 / (M_PI * 6378245.0);
	double tempDouble2 = X * tempDouble1 + Lat0;
	double tempDouble3 = 1 / (cos(tempDouble2 * M_PI / 180));
	dMatData[0][0] = Y * tempDouble1 * tempDouble3 + Lon0;
	dMatData[0][1] = tempDouble2;
}

void CMatrix::GeoToCar(double Lon0, double Lat0, double Lon1, double Lat1)
{
	SetSize(1, 2);
	double tempDouble = M_PI * 6378245.0 / 180;
	dMatData[0][0] = (Lat1 - Lat0) * tempDouble;
	dMatData[0][1] = (Lon1 - Lon0) * tempDouble * cos(Lat1 * M_PI / 180);
}