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
	double **dMatData;    //�������Ԫ�����ݵĶ�ά����
	int iRow;             //�������
	int iCol;             //�������
	// ��������
public:
	int Row() const { return iRow; };                  //������
	int Col() const { return iCol; };                  //������
	void SetSize(int rowN, int colN);                  //��������Ĵ�С��ԭ�����ݲ���
	double& operator () (int row, int col);          //��ȡ����Ԫ��
	double operator () (int row, int col) const;     //���ػ�ȡ����Ԫ�غ�����ֻ��const�����ܷ���
	CMatrix& operator = (const CMatrix& m);
	friend CMatrix operator + (const CMatrix& m1, const CMatrix& m2);    //ע�⣺��Ԫ�������������Լ��ĳ�Ա����
	friend CMatrix operator - (const CMatrix& m1, const CMatrix& m2);
	friend CMatrix operator * (const CMatrix& m1, const CMatrix& m2);
	friend CMatrix operator * (const double& num, const CMatrix& m1);
	friend CMatrix operator * (const CMatrix& m1, const double& num);
	friend CMatrix operator ~ (const CMatrix& m);    //����ת��
	void Unit();      //���ɵ�λ����
	// �������
	void RowExtraction(CMatrix m, int row);       // ����ȡ
	void ColExtraction(CMatrix m, int col);       // ����ȡ
	int* RowOrder(int row);                       // ���ɴ�С����
	int* ColOrder(int col);                       // ���ɴ�С����
	int* MatrixMax();                             // �������ֵ
	// ��ά�ռ�����
	void Vector2D(double X, double Y);                         // ���ɶ�ά������
	void NormVector2D();                                       // ��׼���Ķ�ά������
	void NormVector3D();                                       // ��׼������ά������
	double Vector2DValue();                                    // ��ά����������
	double Vector3DValue();                                    // ��ά����������
	void ACRotate2D(double Phi);                               // ��ά��ת����
	void CRotate2D(double Phi);                                // ��ά��ת����
	double Clock2D(const CMatrix& m1);                         // ������ά��������˳ʱ��нǣ���������붼����Ϊ��ά������
	double AntiClock2D(const CMatrix& m1);                     // ������ά����������ʱ��нǣ���������붼����Ϊ��ά������
	void RegulateVector2(double Angle, double Radius);         // ������ά������������������ͳ���
	// �����滮����
	void Dubins(CMatrix Start, double SD, CMatrix End, double ED, double Rmin);
	void FormationDubins(CMatrix Start, double SD, CMatrix End, double ED, double Rmin, double UAVInitialFormationY, double FormationYLower, double FormationYUpper);
	// �غ����ƶ�����
	void PathMove(double Curvature, double Power);
	// ����ת��
	void CarToGeo(double Lon0, double Lat0, double X, double Y);
	void GeoToCar(double Lon0, double Lat0, double Lon, double Lat);
};

