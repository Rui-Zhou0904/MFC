
// FlightLinePointDlg.h : 头文件
//
#include "Matrix.h"

#pragma once


// CFlightLinePointDlg 对话框
class CFlightLinePointDlg : public CDialogEx
{
// 构造
public:
	CFlightLinePointDlg(CWnd* pParent = NULL);	// 标准构造函数

// 对话框数据
	enum { IDD = IDD_FLIGHTLINEPOINT_DIALOG };

	protected:
	virtual void DoDataExchange(CDataExchange* pDX);	// DDX/DDV 支持


// 实现
protected:
	HICON m_hIcon;

	// 生成的消息映射函数
	virtual BOOL OnInitDialog();
	afx_msg void OnSysCommand(UINT nID, LPARAM lParam);
	afx_msg void OnPaint();
	afx_msg HCURSOR OnQueryDragIcon();
	DECLARE_MESSAGE_MAP()
public:
	double m_EOrigin;
	double m_NOrigin;
	double m_First_EOrigin;
	double m_First_NOrigin;
	double m_Second_EOrigin;
	double m_Second_NOrigin;
	double m_Distance;
	afx_msg void OnBnClickedCal();

	double m_Arrow;
	double m_EOrigin_2;
	double m_NOrigin_2;
	double m_First_Lon;
	double m_First_Lat;
	double m_Distance_2;
	double m_Arrow_2;
	double m_Second_Lon;
	double m_Second_Lat;
	afx_msg void OnBnClickedCal2();
	double m_End_Lon;
	double m_End_Lat;
};
