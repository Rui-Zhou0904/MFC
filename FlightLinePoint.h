
// FlightLinePoint.h : PROJECT_NAME Ӧ�ó������ͷ�ļ�
//

#pragma once

#ifndef __AFXWIN_H__
	#error "�ڰ������ļ�֮ǰ������stdafx.h�������� PCH �ļ�"
#endif

#include "resource.h"		// ������


// CFlightLinePointApp: 
// �йش����ʵ�֣������ FlightLinePoint.cpp
//

class CFlightLinePointApp : public CWinApp
{
public:
	CFlightLinePointApp();

// ��д
public:
	virtual BOOL InitInstance();

// ʵ��

	DECLARE_MESSAGE_MAP()
};

extern CFlightLinePointApp theApp;