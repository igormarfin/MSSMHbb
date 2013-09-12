#ifndef TH1FSyst_h
#define TH1FSyst_h

#include "TH1F.h"
#include "TH2F.h"
#include <iostream>

class TH1FSyst : public TH1F
{
public:
	TH1FSyst():TH1F() { histPlus=0; histMinus=0; histMod=0; };
	TH1FSyst(const TH1F& h2f):TH1F(h2f){};
	~TH1FSyst() {};
	TH1FSyst(TH1F * _histPlus,TH1F * _histMinus,TH1F * _histCentral):TH1F(*_histCentral) { histPlus=_histPlus; histMinus=_histMinus; histCentral=_histCentral;histMod=0; }
	TH1F &  operator()(double _syst) { 




	*histCentral = *this; // restore original vals


/**
temp(pN) =
temp(0) + 1/2(templPlus-temp(0))*pN if pN>0
temp(0) - 1/2(tempMinus-temp(0))*pN if pN<0
**/

		histMod=histCentral;

		if (_syst>0) {

		TH1F _tmpIni( *histPlus);
		_tmpIni.SetName("tmp");
		TH1F * _tmp = new TH1F(_tmpIni);
		_tmp->Add(histCentral,-1);
		_tmp->Scale(1./2.);
		 histMod->Add( _tmp,1*_syst);
		delete _tmp;
		}

		if (_syst<0) { 
		TH1F _tmpIni( *histMinus);
		_tmpIni.SetName("tmp");
                TH1F * _tmp = new TH1F(_tmpIni);
                _tmp->Add(histCentral,-1);
                _tmp->Scale(1./2.);
                 histMod->Add( _tmp,-1*_syst);
                delete _tmp;
                }

	return *histMod;

	 }
private:

TH1F * histPlus;
TH1F *histMinus;
TH1F *histCentral;
TH1F * histMod;

};


class TH2FSyst : public TH2F
{
public:
	TH2FSyst():TH2F() { histPlus=0; histMinus=0; histMod=0; };
	TH2FSyst(const TH2F& h2f):TH2F(h2f){};
	~TH2FSyst() {};
	TH2FSyst(TH2F * _histPlus,TH2F * _histMinus,TH2F * _histCentral):TH2F(*_histCentral) { histPlus=_histPlus; histMinus=_histMinus; histCentral=_histCentral;histMod=0; }
	TH2F &  operator()(double _syst) { 




	*histCentral = *this; // restore original vals


/**
temp(pN) =
temp(0) + 1/2(templPlus-temp(0))*pN if pN>0
temp(0) - 1/2(tempMinus-temp(0))*pN if pN<0
**/

		histMod=histCentral;

		if (_syst>0) {

		TH2F _tmpIni( *histPlus);
		_tmpIni.SetName("tmp");
		TH2F * _tmp = new TH2F(_tmpIni);
		_tmp->Add(histCentral,-1);
		_tmp->Scale(1./2.);
		 histMod->Add( _tmp,1*_syst);
		delete _tmp;
		}

		if (_syst<0) { 
		TH2F _tmpIni( *histMinus);
		_tmpIni.SetName("tmp");
                TH2F * _tmp = new TH2F(_tmpIni);
                _tmp->Add(histCentral,-1);
                _tmp->Scale(1./2.);
                 histMod->Add( _tmp,-1*_syst);
                delete _tmp;
                }

	return *histMod;

	 }
private:

TH2F * histPlus;
TH2F *histMinus;
TH2F *histCentral;
TH2F * histMod;

};




#endif

