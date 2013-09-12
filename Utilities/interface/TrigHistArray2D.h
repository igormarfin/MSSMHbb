#ifndef TrigHistArray2D_h
#define TrigHistArray2D_h


class TrigHistArray2D {
public:
  std::vector<TH2F*> histos;
  TH2F* histAllTrig;
  TH2F* histAllTrigWeighted;
  std::vector<std::string>* gtl;  // genericTriggerList
  unsigned int filterPattern;
  
  TrigHistArray2D() {
    std::cout << "Standard constructor called, this should not happen" << std::endl;
  }

  TrigHistArray2D(std::vector<std::string>* theGtl,std::vector<std::string>* theFilter,const char* genName,const char* genTitle,int nx,float xmin,float xmax,int ny,float ymin,float ymax) {
    gtl = theGtl;
    char ctn[1000];
    char ctt[1000];
    // first set the filter pattern
    filterPattern = 0;
    if (theFilter != NULL) {
      for (std::vector<std::string>::iterator fit=theFilter->begin(); fit != theFilter->end(); ++fit) {
	std::vector<std::string>::iterator tSlotBtag = std::find(theGtl->begin(), theGtl->end(), *fit);
	if (tSlotBtag != theGtl->end()) {
	  filterPattern = filterPattern | (1<<(tSlotBtag - theGtl->begin()));
	} else {
	  std::cout << "Filter trigger " << *fit << " not found in any slot" << std::endl;
	}
      }
    } else {
      filterPattern = ~0;  // set all filter bits to one
    }

    for (unsigned int ib=0; ib<gtl->size(); ++ib) {
      TH2F* theHist = 0;
      if (filterPattern & (1<<ib)) {
	sprintf(ctn,"%sTrig%d",genName,ib);
	sprintf(ctt,"%s Trigger %s",genTitle,(*gtl)[ib].c_str());
	theHist = new TH2F(ctn,ctt,nx,xmin,xmax,ny,ymin,ymax);
      }
      histos.push_back(theHist);
    }
    // inclusive histogram
    sprintf(ctn,"%sAllTrig",genName);
    sprintf(ctt,"%s all triggers",genTitle);
    histAllTrig = new TH2F(ctn,ctt,nx,xmin,xmax,ny,ymin,ymax);
    sprintf(ctn,"%sAllTrigWeighted",genName);
    sprintf(ctt,"%s all triggers weighted",genTitle);
    histAllTrigWeighted = new TH2F(ctn,ctt,nx,xmin,xmax,ny,ymin,ymax);
  }
			     
  void fill(int theTrgAccept,float x,float y,float weight=1) {
    for (unsigned int ib=0; ib<gtl->size(); ++ib) {
      if (theTrgAccept & filterPattern & (1<<ib)) {
	if (histos[ib] != NULL) {
	  histos[ib]->Fill(x,y,weight);
	} else {
	  std::cout << "TrigHistArray2D::Fill: bad histogram pointer " << std::endl;
	}
      }
    }
    histAllTrig->Fill(x,y);
    histAllTrigWeighted->Fill(x,y,weight);
  }

  void fillTW(float x,float y,float* trigWeight,float weight=1) {
    // instead of trigger selection, fill with external trigger efficiency weight
    for (unsigned int ib=0; ib<gtl->size(); ++ib) {
      if (filterPattern & (1<<ib)) {
	if (histos[ib] != NULL) {
	  histos[ib]->Fill(x,y,trigWeight[ib]*weight);
	} else {
	  std::cout << "TrigHistArray2D::Fill: bad histogram pointer " << std::endl;
	}
      }
    }
    histAllTrig->Fill(x,y);
    histAllTrigWeighted->Fill(x,y,weight);
  }

  static void mergeHistos( const int nMerge, TrigHistArray2D* taMerge[] ) {
    // merge histograms [0]... [nMerge-1], store into [nMerge]
    if ( taMerge[nMerge] == NULL ) {
      std::cout << "TrigHistArray::mergeHistos: Cannot find histograms " << std::endl;
      return;
    }
    for (unsigned int ihist=0; ihist<taMerge[nMerge]->histos.size(); ++ihist) {
      if (taMerge[nMerge]->histos[ihist] != NULL) {
	for (int itpat=0; itpat<3; ++itpat) {
	  if (taMerge[itpat]->histos[ihist] != NULL) {
	    taMerge[nMerge]->histos[ihist]->Add( taMerge[itpat]->histos[ihist] );
	  } else {
	    std::cout << "TrigHistArray::mergeHistos: Cannot find hist for itpat=" << itpat << std::endl;
	  }
	}
	// replace square-summed errors by linear sum to account for correlation
	for (int ibinx=1; ibinx<=taMerge[nMerge]->histos[ihist]->GetXaxis()->GetNbins(); ++ibinx) {
	  for (int ibiny=1; ibiny<=taMerge[nMerge]->histos[ihist]->GetYaxis()->GetNbins(); ++ibiny) {
	    float theError = 0;
	    for (int itpat=0; itpat<3; ++itpat) {
	      theError += taMerge[itpat]->histos[ihist]->GetBinError( ibinx,ibiny );
	    }
	    taMerge[nMerge]->histos[ihist]->SetBinError( ibinx, ibiny, theError );
	  }
	}
      }
    }
  }

  void add(const TrigHistArray2D* a1, const TrigHistArray2D* a2, double c1 = 1, double c2 = 1) {
    std::cout << "TrigHistArray2D: a1=" << a1 << "  a2=" << a2 << std::endl;
    for (unsigned int ib=0; ib<gtl->size(); ++ib) {
      if (filterPattern & (1<<ib)) {
	if ( (a1->histos[ib] != NULL) && (a2->histos[ib] != NULL)) {
	  histos[ib]->Add( a1->histos[ib], a2->histos[ib], c1, c2 );
	} else {
	  std::cout << "TrigHistArray2D::Fill: bad histogram pointer " << std::endl;
	}
      }
    }
    histAllTrig->Add( a1->histAllTrig, a2->histAllTrig, c1, c2 );
    histAllTrigWeighted->Add( a1->histAllTrigWeighted, a2->histAllTrigWeighted, c1, c2 );
    return;
  }
};

#endif // #ifdef 

