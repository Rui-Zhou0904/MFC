
// FlightLinePointDlg.h : ͷ�ļ�
//
#include "Matrix.h"

#pragma once


// CFlightLinePointDlg �Ի���
class CFlightLinePointDlg : public CDialogEx
{
// ����
public:
	CFlightLinePointDlg(CWnd* pParent = NULL);	// ��׼���캯��

// �Ի�������
	enum { IDD = IDD_FLIGHTLINEPOINT_DIALOG };

	protected:
	virtual void DoDataExchange(CDataExchange* pDX);	// DDX/DDV ֧��


// ʵ��
protected:
	HICON m_hIcon;

	// ���ɵ���Ϣӳ�亯��
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
