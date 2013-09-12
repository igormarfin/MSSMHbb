#include "lsPrint.h"

#include <iostream>

using namespace std;

void PrintHelp(){
  cout << "lsRoot \n"; 
  cout << "Author Justin Griffiths griffith@cern.ch\n";
  cout << "With thanks to the root development team\n";
  cout << "lsRoot Basic functionality:\n"; 
  cout << "    lsRoot file Print all TObjects in baseDir of file (File need not end .root!)\n";
  cout << "    lsRoot dir/ (run program on all root files in directory) \n";
  cout << "    lsRoot -R file.root (print out root file recursively) \n";
  cout << "    lsRoot -Z dir/ (run program on all root files in directory recursively) \n";
  cout << "    lsRoot -p */Muon/MuonR*/MDT/* file.root (start printout in rootPath '*/Muon/MuonR*/MDT/*') \n";
  cout << "    lsRoot -p regex '.*/Muon/Muon[RS].*/[MTB].*' file.root (start printout in rootPath '.*/Muon/Muon[RS].*/[MTB].*') \n";
  cout << "    lsRoot -l file.root (print in long format)\n";
  cout << "    lsRoot -ltrh file.root (print in long format sort by reverse time nBytes in human readable format)\n";
  cout << "    lsRoot -C file.root (attempts to approximate bare 'ls' functionality on by default)\n";
  cout << "    lsRoot -S ouput h_nTrk -p '*/lb_*/Tau/*' file.root (create output.ps with all hists w/ name h_nTrk encountered\n";
  cout << "            can be list of hists 'lsRoot --help script')\n";
  cout << "    lsRoot --draw  Draw TH1/TH2s in root file (control which to draw w/ -p path -R (recurseive or not))\n\n";

  cout << "General Usage:\n";
  cout << "lsRoot should be able to receive and give to pipe/redirects i.e. 'ls -d dir/*pattern* | lsRoot -args | grep \"Keys\" \n";
  cout << " -A              Print Mean of TH1\n";
  cout << " --alpha         Sort keys alphabetically (default)\n";
  cout << " -B              Print Nbins of TH1\n";
  cout << " -c/--custom     For Each TKey choose your own layout '\\h,\\I,\\B,\\A,\\N,\\D,\\d,\\c' more w/ 'lsRoot --help custom' \n";
  cout << " --cycle [opt]   opt = 'all' print all cycles of like keys, 'last' print last cycle of all keys, 'n' print nth cycle of keys.  Default 'last'\n";
  cout << " -C              Display in a manner similar to ls (multiColumn) default if to terminal\n";
  cout << " -d              Print TDirectory path w/ each key\n";
  cout << " -dd             Print bashPath/rootfile:/TDirectory path w each key\n";
  cout << " --draw          Draw Ascii version of histogram (TH1 or TH2) only 'lsRoot --help draw' for more info\n";
  cout << " -D              Print RMS of TH1\n";
  cout << " -e              Show root Warning messages Default off\n";
  cout << " -h              human readable\n";
  cout << " --help          Print this help message\n";
  cout << " -i              same as 'lsRoot --draw --interactive' 'lsRoot --help draw' for more info\n"; 
  cout << " --interactive   assumes '--draw' the program is paused after each histogram is drawn where the user can then redraw the histogram changing axis/log mode\n";
  cout << " -I              Print TH1::Integral() w/ each TH1\n";
  cout << " -l              List in long format\n";
  cout << " -m              Map the file \n";
  cout << " -N              Print Number of entries of TH1 do '-N[bin]' to print only content in certain bin\n";
  cout << " --noColor/ -G   Do not display color (color may produce perverse results in some terminals) Disabled in Windows\n";
  cout << " -p rootPath     Give intial root path to start looping through *,? accepted, TPRegex if -p regex 'lsRoot --help path'\n";
  cout << " -P rootPath     Same as 'lsRoot -p regex'\n";
  cout << " --precision v   Set the precision of -IBAND or --custom '\\I\\B\\A\\N\\D' outputs\n";
  cout << " -r              Reverse order keys\n";
  cout << " -R              Loop through individual root files recursively\n";
  cout << " -S fName hList  'lRoot --help script' for more info \n";
  cout << " --size          Sort TObjs by size\n";
  cout << " --sortN         Sort TObjs by (TH1)->GetEntries (0 if no cast works i.e. TDirectory is 0)\n";
  cout << " --sortNBins     Sort TObjs by (TH1)->GetNbinX()*Y()*Z() (0 if no cast works i.e. TDirectory is 0)\n";
  cout << " --stats         count recursively nKeys, nBytes, nBins for each TDirectory (combine w/ -R to get)\n";
  cout << " -t              Sort keys by date (like ls -t) is faster (even faster -rt)! (default is alphabetical)\n";
  cout << " -T              Print Any and all trees w/ opt topOnly\n";
  cout << " -TT             Print Any and all trees w/o opt topOnly\n";
  cout << " -TTT            Print Any and all trees recursively w/ custom lsRoot verbose mode includes type of variable\n";
  cout << " --tree 'code'   Print code to stdout where '\\n' is new line '\\t' is branch type '\\b' is branch name\n";
  cout << " --types 'list'  comma separated list of TClasses such as 'TH1F,TH2F,..'\n";
  cout << " -Z              Loop through all root files in bash directory recursively \n";
  cout << "\nTo use color formatting 'cp .lsRoot_colors to ~/' see file to set user defined colors see '--help color' \n\n";
}

void PrintScriptHelp(){
  cout << "Option of lsRoot to collect all objects castable to TH1 (TH[123][A-Z] TProfile[123][A-Z] etc.) and print them to\n";
  cout << "a post-script file.  Also dump directly TCanvases(containing TH1/THStack/TGraph castable objects) into ps.\n";
  cout << "Basic Usage:\n";
  cout << "  lsRoot -S myOutFile 'hist1,hist2,...' file.root\n";
  cout << "  The above will create myOutFile.ps which will draw per page hist1,hist2,... up to 8 histograms per page. \n";
  cout << "  A user can then call pstopdf(mac) or ps2pdf(linux) to convert myOutFile.ps to pdf.\n";
  cout << "  If the user wants to directly create the pdf name 'myOutFile'-->'myOutFile.pdf'\n";
  cout << "  lsRoot will then determine if first command pstopdf exists on your machine then ps2pdf\n";
  cout << "  If one of the commands is found, then lsRoot will convert the ps to a pdf and remove the post-script file\n\n";

  cout << "  The hist list is a comma separated list where histn can be pure strings or regular expressions\n";
  cout << "  understood by TPRegexp (i.e. .*,[A-Z]*,[0-5]+,\\w,(str1|str2), etc.).  If you want histn to be drawn in\n";
  cout << "  logy/z write histName'logy/z' see examples below for precise syntax(only for pdf).  The idea is if you have a root\n";
  cout << "  file with a rich directory structure such as path1/path2/[dummy]/path4 where in each of path 4 there are\n";
  cout << "  histograms you'd like to see altogether and quickly you can combine '-S fName hList' w/ \n";
  cout << "  '-p */dir1/dir2/*/dir4/ to gather these hists into a pdf. \n\n";

  cout << " Example usage:\n";
  cout << "  lsRoot -p */Muon/*/MDT/MDT*/Chambers/BIL* -S BILXXXXChamberStudy '.*MDT_Station_TUBE_ADCCutLogy,.*MDT_Station_TDCADC.*Logz,.*MDT_Station_LAY.*' Monitor.root\n";
  cout << "  lsRoot -p */Muon/*/MDT/Overview -S OverviewDump '.*'  Monitor.root (This dumps all hists found in directory to pdf w/ one per page)\n";
  cout << "  \n";
}

void PrintPathHelp(){
  cout << "Option to start looping through root file at specified path *,? accepted\n";
  cout << "Basic usage:\n";
  cout << "  lsRoot -p path file.root (looks for path 'path' in root file then starts list from there)\n";
  cout << "  lsRoot -p dir1/*/dir2/*ir3/d*r4/* file.root \n";
  cout << "\nNotice Difference Below\n";
  cout << "  lsRoot -p dir1/so?eth*g/*partialHistName* file.root (The last portion will be considered TObject (includes TDirectory))\n";
  cout << "  lsRoot -p dir1/so?eth*g/*partialHistName*/ file.root (The last portion will be considered TDirectory only)\n";
  cout << "  lsRoot -p regex 'anything accepted by TPRegexp' \n";
}

void PrintDrawHelp(){
  cout << "Draw Ascii versions of TH1 and TH2 hists\n";
  cout << "Output graphics are configurable in your ~/.lsRoot_colors.\n";
  cout << "Note that the window size in Windows is assumed to be 80Rowsx60Cols\n";
  cout << "This can be altered in .lsRoot_colors 'nRows xx' 'nCols xx' \n\n";
  
  cout << "    StatBox color\n";
  cout << "    XAxis color\n";
  cout << "    YAxis color\n";
  cout << "    FillColor color\n";
  cout << "    Here color refers to any of the escape codes as provided by running 'lsRoot --help color'\n";
  cout << "    SetPalette color1, color2, color3, ...\n";
  cout << "    Here you can define a color palette for TH2 histograms by listing any sized list of colors as displayed from 'lsRoot --help color'\n\n";
  cout << "    FillSymbol X0\n";
  cout << "    a string of length 1 or 2(if longer only first two characters are used) this is the ascii character used\n";
  cout << "    as the FillStyle of the histogram, when length is 2, then all adjacent bins will have a different 'FillSymbol'\n";
  cout << "    This is useful to distinguish b/w differing bins esp. if bins have nearly the same height\n\n";
  cout << "    SetOptStat 1111 (Not yet implemented)\n\n";
  cout << "    LogyThreshold 800 (if h->GetMaximum/h->GetMinimum(0.000001) > LogyThreshold) then TH1 is in logMode\n";
  cout << "    TH2 Logz is applied if h->GetMaximum/h->GetMaximum(0.000001) > 100 (not configurable)\n\n";   
  cout << "    if 'lsRoot -i' or 'lsRoot --draw --interactive' or 'lsRoot --interactive' Draw is interactive\n";
  cout << "    The histogram is drawn as normal after which the program pauses.\n";
  cout << "    The program will resume w/ a 'return' command, or see more options w/ 'i' to eneter the interactive mode.\n";
  cout << "    At this point one can enter 'x','y', or 'z' to change the axis range.  'l' to toggle b/w log mode. 'd' to redraw the hist\n";
  cout << "    'r' to reset all previous changes, or 'q' to continue to the next histogram.\n";
  cout << "    All changes are cummulative in interactive mode.  That is if you change the x-axis range and then draw.  Then change the y-axis and then\n";
  cout << "    redraw again, the resultant hist will have both the x and y axises resized.\n\n";
}

void PrintCustomHelp(){
  cout << "Draw each TKey w/ a user defined custom format\n";
  cout << "Run w/ 'lsRoot -c/--custom string'\n";
  cout << "i.e. to print keyName w/ class type: \n";
  cout << "  lsRoot --custom '\\h \\c'\n\n";
  cout << "The available escape characters:\n";
  cout << "   \\n              New Line\n";
  cout << "   \\p              strip pad formating from specific line\n";
  cout << "   \\h              TKey::GetName()\n";
  cout << "   \\c              TKey::GetClassName()\n";
  cout << "   \\C              TKey::GetCycle()\n";
  cout << "   \\d              TDirectory::GetPath() (bash portion stripped)\n";
  cout << "   \\N              TH1/TTree::GetEntries()\n";
  cout << "   \\N[bin]         TH1::GetBinContent(bin)\n";
  cout << "   \\N[binx-biny]   TH1::GetBinContent(bin) on each bin in range\n";
  cout << "   \\N[all]         TH1::GetBinContent(bin) on all bins\n";
  cout << "   \\I              TH1::Integral()\n";
  cout << "   \\I[x1-x2]       TH1::Integral(x1,x2)\n";
  cout << "   \\I[x1-x2,y1-y2] TH2::Integral(x1,x2,y1,y2)\n";
  cout << "   \\A              TH1::GetMean()\n";
  cout << "   \\AX             TH1::GetMean(1)\n";
  cout << "   \\AY             TH1::GetMean(2)\n";
  cout << "   \\AZ             TH1::GetMean(3)\n";
  cout << "   \\D              TH1::GetRMS()\n";
  cout << "   \\DX             TH1::GetRMS(1)\n";
  cout << "   \\DY             TH1::GetRMS(2)\n";
  cout << "   \\DZ             TH1::GetRMS(3)\n";
  cout << "   \\B              TH1::GetNbinsX()*TH1::GetNbinsY()*TH1::GetNbinsZ()\n";
  cout << "   \\BX             TH1::GetNbinxX()\n";
  cout << "   \\BY             TH1::GetNbinxY()\n";
  cout << "   \\BZ             TH1::GetNbinxZ()\n";
  cout << "   \\t              TKey::GetDaTime()->AsString()\n";
  cout << "   \\s              TKey::GetNbytes()\n";
  cout << "   \\k              TKey::GetNkeys()\n";
  cout << "   expr<TFormula>  Evaluate an expression such as TH1::GetEntries()*TH1::GetMean() 'expr<\\N*\\A>'\n";
  cout << "                   expressions are evaluated at the end such that all \\X where X is an escape\n";
  cout << "                   charcter described above.  For example if lsRoot --custom 'string' produces\n";
  cout << "                   output that appears to be evaluable, then lsRoot --custom 'expr<string>' \n";
  cout << "                   Will evaluate it.\n";
  cout << "                   One can save their favorite custom formats in the ~/.lsRoot_colors file such as:\n\n";
  cout << "                   customFormat myFormat =\\A\\B\\N\\D\\I \\d/\\h \n\n";
  cout << "                   Then run simply run w/ 'lsRoot -c myFormat file[.root]'\n";
  cout << endl;
  cout << "One can use this feature to generate simple code such as this:\n";
  cout << "lsRoot --custom '\\p\\c* \\h = (\\c*) _file0->Get(\"\\d\\h\");' file.root:/*/dir/*/hist*\n";
  cout << endl;
}

void showColors(){

  cout << "To Use Color simply place .lsRoot_colors in $HOME directory.  \n";
  cout << "In .lsRoot_colors to make all TDirectories blue do:\n";
  cout << "TDirectoryFile 01;34\n\n";    
  
  cout << "\033[0;00m" << "00;00\n";
  cout << "\033[0;30m" << "00;30\n";
  cout << "\033[1;30m" << "01;30\n";
  cout << "\033[4;30m" << "04;30\n";
  cout << "\033[5;30m" << "05;30\n";
  cout << "\033[0;31m" << "00;31\n";
  cout << "\033[1;31m" << "01;31\n";
  cout << "\033[4;31m" << "04;31\n";
  cout << "\033[5;31m" << "05;31\n";
  cout << "\033[0;32m" << "00;32\n";
  cout << "\033[1;32m" << "01;32\n";
  cout << "\033[4;32m" << "04;32\n";
  cout << "\033[5;32m" << "05;32\n";
  cout << "\033[0;33m" << "00;33\n";
  cout << "\033[1;33m" << "01;33\n";
  cout << "\033[4;33m" << "04;33\n";
  cout << "\033[5;33m" << "05;33\n";
  cout << "\033[0;34m" << "00;34\n";
  cout << "\033[1;34m" << "01;34\n";
  cout << "\033[4;34m" << "04;34\n";
  cout << "\033[5;34m" << "05;34\n";
  cout << "\033[0;35m" << "00;35\n";
  cout << "\033[1;35m" << "01;35\n";
  cout << "\033[4;35m" << "04;35\n";
  cout << "\033[5;35m" << "05;35\n";
  cout << "\033[0;36m" << "00;36\n";
  cout << "\033[1;36m" << "01;36\n";
  cout << "\033[4;36m" << "04;36\n";
  cout << "\033[5;36m" << "05;36\n";
  cout << "\033[0;37m" << "00;37\n";
  cout << "\033[1;37m" << "01;37\n";
  cout << "\033[4;37m" << "04;37\n";
  cout << "\033[5;37m" << "05;37\n";
  cout << "\033[0;38m" << "00;38\n";
  cout << "\033[1;38m" << "01;38\n";
  cout << "\033[4;38m" << "04;38\n";
  cout << "\033[5;38m" << "05;38\n";
  cout << "\033[0;39m" << "00;39\n";
  cout << "\033[1;39m" << "01;39\n";
  cout << "\033[4;39m" << "04;39\n";
  cout << "\033[5;39m" << "05;39\n";
}


