
// FlightLinePointDlg.cpp : 实现文件
//

#include "stdafx.h"
#include "FlightLinePoint.h"
#include "FlightLinePointDlg.h"
#include "afxdialogex.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#endif
#define M_PI 3.1415926

// 用于应用程序“关于”菜单项的 CAboutDlg 对话框

class CAboutDlg : public CDialogEx
{
public:
	CAboutDlg();

// 对话框数据
	enum { IDD = IDD_ABOUTBOX };

	protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV 支持

// 实现
protected:
	DECLARE_MESSAGE_MAP()
};

CAboutDlg::CAboutDlg() : CDialogEx(CAboutDlg::IDD)
{
}

void CAboutDlg::DoDataExchange(CDataExchange* pDX)
{
	CDialogEx::DoDataExchange(pDX);
}

BEGIN_MESSAGE_MAP(CAboutDlg, CDialogEx)
END_MESSAGE_MAP()


// CFlightLinePointDlg 对话框



CFlightLinePointDlg::CFlightLinePointDlg(CWnd* pParent /*=NULL*/)
	: CDialogEx(CFlightLinePointDlg::IDD, pParent)
	, m_EOrigin(105)
	, m_NOrigin(39)
	, m_First_EOrigin(105)
	, m_First_NOrigin(39)
	, m_Second_EOrigin(105)
	, m_Second_NOrigin(39)
	, m_Distance(0)
	, m_Arrow(0)
	, m_EOrigin_2(105)
	, m_NOrigin_2(39)
	, m_First_Lon(105)
	, m_First_Lat(39)
	, m_Distance_2(0)
	, m_Arrow_2(0)
	, m_Second_Lon(105)
	, m_Second_Lat(40)
	, m_End_Lon(0)
	, m_End_Lat(0)
{
	m_hIcon = AfxGetApp()->LoadIcon(IDR_MAINFRAME);
}

void CFlightLinePointDlg::DoDataExchange(CDataExchange* pDX)
{
	CDialogEx::DoDataExchange(pDX);
	DDX_Text(pDX, IDC_EDIT6, m_EOrigin);
	DDX_Text(pDX, IDC_EDIT7, m_NOrigin);
	DDX_Text(pDX, IDC_EDIT1, m_First_EOrigin);
	DDX_Text(pDX, IDC_EDIT2, m_First_NOrigin);
	DDX_Text(pDX, IDC_EDIT3, m_Second_EOrigin);
	DDX_Text(pDX, IDC_EDIT4, m_Second_NOrigin);
	DDX_Text(pDX, IDC_EDIT5, m_Distance);
	DDX_Text(pDX, IDC_EDIT8, m_Arrow);
	DDX_Text(pDX, IDC_EDIT14, m_EOrigin_2);
	DDX_Text(pDX, IDC_EDIT15, m_NOrigin_2);
	DDX_Text(pDX, IDC_EDIT9, m_First_Lon);
	DDX_Text(pDX, IDC_EDIT10, m_First_Lat);
	DDX_Text(pDX, IDC_EDIT13, m_Distance_2);
	DDX_Text(pDX, IDC_EDIT16, m_Arrow_2);
	DDX_Text(pDX, IDC_EDIT11, m_Second_Lon);
	DDX_Text(pDX, IDC_EDIT12, m_Second_Lat);
	DDX_Text(pDX, IDC_EDIT17, m_End_Lon);
	DDX_Text(pDX, IDC_EDIT18, m_End_Lat);
}

BEGIN_MESSAGE_MAP(CFlightLinePointDlg, CDialogEx)
	ON_WM_SYSCOMMAND()
	ON_WM_PAINT()
	ON_WM_QUERYDRAGICON()
	ON_BN_CLICKED(IDC_CAL, &CFlightLinePointDlg::OnBnClickedCal)
	ON_BN_CLICKED(IDC_CAL2, &CFlightLinePointDlg::OnBnClickedCal2)
END_MESSAGE_MAP()


// CFlightLinePointDlg 消息处理程序

BOOL CFlightLinePointDlg::OnInitDialog()
{
	CDialogEx::OnInitDialog();

	// 将“关于...”菜单项添加到系统菜单中。

	// IDM_ABOUTBOX 必须在系统命令范围内。
	ASSERT((IDM_ABOUTBOX & 0xFFF0) == IDM_ABOUTBOX);
	ASSERT(IDM_ABOUTBOX < 0xF000);

	CMenu* pSysMenu = GetSystemMenu(FALSE);
	if (pSysMenu != NULL)
	{
		BOOL bNameValid;
		CString strAboutMenu;
		bNameValid = strAboutMenu.LoadString(IDS_ABOUTBOX);
		ASSERT(bNameValid);
		if (!strAboutMenu.IsEmpty())
		{
			pSysMenu->AppendMenu(MF_SEPARATOR);
			pSysMenu->AppendMenu(MF_STRING, IDM_ABOUTBOX, strAboutMenu);
		}
	}

	// 设置此对话框的图标。  当应用程序主窗口不是对话框时，框架将自动
	//  执行此操作
	SetIcon(m_hIcon, TRUE);			// 设置大图标
	SetIcon(m_hIcon, FALSE);		// 设置小图标

	// TODO:  在此添加额外的初始化代码

	return TRUE;  // 除非将焦点设置到控件，否则返回 TRUE
}

void CFlightLinePointDlg::OnSysCommand(UINT nID, LPARAM lParam)
{
	if ((nID & 0xFFF0) == IDM_ABOUTBOX)
	{
		CAboutDlg dlgAbout;
		dlgAbout.DoModal();
	}
	else
	{
		CDialogEx::OnSysCommand(nID, lParam);
	}
}

// 如果向对话框添加最小化按钮，则需要下面的代码
//  来绘制该图标。  对于使用文档/视图模型的 MFC 应用程序，
//  这将由框架自动完成。

void CFlightLinePointDlg::OnPaint()
{
	if (IsIconic())
	{
		CPaintDC dc(this); // 用于绘制的设备上下文

		SendMessage(WM_ICONERASEBKGND, reinterpret_cast<WPARAM>(dc.GetSafeHdc()), 0);

		// 使图标在工作区矩形中居中
		int cxIcon = GetSystemMetrics(SM_CXICON);
		int cyIcon = GetSystemMetrics(SM_CYICON);
		CRect rect;
		GetClientRect(&rect);
		int x = (rect.Width() - cxIcon + 1) / 2;
		int y = (rect.Height() - cyIcon + 1) / 2;

		// 绘制图标
		dc.DrawIcon(x, y, m_hIcon);
	}
	else
	{
		CDialogEx::OnPaint();
	}
}

//当用户拖动最小化窗口时系统调用此函数取得光标
//显示。
HCURSOR CFlightLinePointDlg::OnQueryDragIcon()
{
	return static_cast<HCURSOR>(m_hIcon);
}



void CFlightLinePointDlg::OnBnClickedCal()
{
	// TODO:  在此添加控件通知处理程序代码
	UpdateData(1);
	double EOrigin = m_EOrigin;
	double NOrigin = m_NOrigin;

	double TargetAreaLon = m_First_EOrigin;
	double TargetAreaLat = m_First_NOrigin;
	double dTargetLon = m_Second_EOrigin;
	double dTargetLat = m_Second_NOrigin;

	CMatrix PointO(1, 2); CMatrix PointA(1, 2); CMatrix PointAO(1, 2);
	PointO.GeoToCar(EOrigin, NOrigin, TargetAreaLon, TargetAreaLat);
	PointA.GeoToCar(EOrigin, NOrigin, dTargetLon, dTargetLat);
	PointAO = PointA - PointO;
	double distanceAO = sqrt(PointAO(0, 0)*PointAO(0, 0) + PointAO(0, 1)*PointAO(0, 1));
	m_Distance = distanceAO;
	double BasicArrowAngle;
	CMatrix tempArrowElement(1, 2);
	PointAO.NormVector2D();
	tempArrowElement(0, 0) = 1.0; tempArrowElement(0, 1) = 0.0;
	BasicArrowAngle = (180 / M_PI) * tempArrowElement.AntiClock2D(PointAO);//弧度
	m_Arrow=BasicArrowAngle;
	UpdateData(0);
}


void CFlightLinePointDlg::OnBnClickedCal2()
{
	// TODO:  在此添加控件通知处理程序代码
	UpdateData(1);
	double EOrigin_2 = m_EOrigin_2;
	double NOrigin_2 = m_NOrigin_2;
	double TargetAreaLon_2 = m_First_Lon;
	double TargetAreaLat_2 = m_First_Lat;
	double dTargetLon_2 = m_Second_Lon;
	double dTargetLat_2 = m_Second_Lat;
	double Distance_2 = m_Distance_2;
	double Arrow_2 = m_Arrow_2;

	CMatrix PointO2(1, 2); CMatrix PointA2(1, 2); CMatrix PointAO2(1, 2);
	PointO2.GeoToCar(EOrigin_2, NOrigin_2, TargetAreaLon_2, TargetAreaLat_2);
	PointA2.GeoToCar(EOrigin_2, NOrigin_2, dTargetLon_2, dTargetLat_2);
	PointAO2 = PointA2 - PointO2;
	CMatrix TempArrow1(PointAO2);
	TempArrow1.RegulateVector2(Arrow_2, Distance_2);
	CMatrix TempArrow2(1, 2);
	double a = TempArrow1(0, 0);
	double b = TempArrow1(0, 1);
	TempArrow2.CarToGeo(EOrigin_2, NOrigin_2, a, b);
	m_End_Lon = TempArrow2(0, 0);
	m_End_Lat = TempArrow2(0, 1);
	UpdateData(0);
}
