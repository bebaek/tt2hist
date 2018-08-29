/*

Make histogram files from time-tagged data files (input files in .ht2 by hydraharp)
BB, 2008

v21:
process SYNC; ch 0 becomes sync.

v23c:
port to vc++ express 2010

*/
#define VERSION "24"		// version of this program

#include "ht2read.h"		// ht2 file handling class
#include <windows.h>
#include <assert.h>
#include <stdio.h>
#include <conio.h>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <ctime>			// measuring code execution time
#include <iomanip>			// c++ output formatting
#include <cmath>				// using pow()
#include <algorithm>		// sort function

//#define BUFFER_SIZE 0xFE00
//#define NTTMAX		0x400000	// timestamp data size per processing
//#define NTTMAX		1000000	// timestamp data size per processing. look for this in ht2read.h
#define TCONV		1e-12		// hardware time bin size
#define I2CH(c)		c			// array index to channel designator mapper
#define CH2I(c)		c			// "internal" channel designator to array index mapper
//#define I2CHINV(c)	c-1
#define NCHTOT		8			// n of channels including copied ones.
using namespace std;

struct datastruct { char channel; long long timestamp; };

// main correlation calculating class
class corr
{
	double w;				// time window
	unsigned long long wbin;
public:
	unsigned int *c2[NCHTOT][NCHTOT];			// count arrays of 2nd order correlation
	//unsigned int *c3;				// count arrays of 3rd order correlation
	unsigned int *hc;				// count arrays of the higher order (>2) correlation

	unsigned long long ndatareq;	// number of data to read at one time.
	char *ftt;				// data file name
	datastruct *data;		// main data
	long long binsize2;			// number of merged bins for 2-fold
	long long hbinsize;			// number of merged bins for 3-fold
	long nbin2;				// number of bins for 0..boundary for 2-fold
	long hnbin;				// number of bins for 0..boundary for 3-fold
	int ch[NCHTOT];				// channels to process for high-order hist
	int highorder;			// the order of correlation to calculate: 3 or 4
	long long ttoffset[NCHTOT];	// timetag offset
	long long ttoffsetmin, ttoffsetmax;
	int wantoffset;			// flag for offset processing
	int stopafterround;		// round # to stop after
	int wantstopafter;		// flag for stop after round #
	char *fpre;
	int rep;				// desired repetition input
	int irep;				// repetition monitor
	unsigned long long ctot[NCHTOT];		// total counts
	unsigned long long ttot;		// measurement time resulted
	unsigned long long cbyround[NCHTOT];		// total ch counts in a round
	unsigned long long tbyround;		// measurement time in a round
	unsigned long long c2totbyround[NCHTOT][NCHTOT];	// total hist2 counts in a round
	unsigned long long hctotbyround;	// total hist3 (or 4) counts in a round
	unsigned int ncarryover;	// number of carried-over data
	unsigned int ndata;
	unsigned int ndataprev;
	unsigned long long tstart0;		// very-start timestamp
	unsigned int istart;		// starting data index at the current data set
	unsigned int iend;			// ending data index at the current data set
	double norm2[NCHTOT][NCHTOT];				// normalization for g2
	double hnorm;
	// things for hhist
	long **tb;
	unsigned int *ndet;
	int noheader;				// indicate no header in ht2 file
	
	// variables introduced in v18: piecewise meas values
	unsigned int *c2byround[NCHTOT][NCHTOT];			// count arrays of 2nd order correlation
	unsigned int *hcbyround;					// count arrays of the higher order (>2) correlation
	double norm2byround[NCHTOT][NCHTOT];				// normalization for g2
	double hnormbyround;
	double *g2byround[NCHTOT][NCHTOT];		// g2 in a round
	double *hgbyround;					// high-order g in a round
	double *g2avg[NCHTOT][NCHTOT];			// g2 average (save this)
	double *hgavg;						// high-order g average (save this)
	unsigned int hnbintot;				// total number of bins in highorder histogram

	// for datacopying
	datastruct *datacopy;				// datacopy
	int wantcopy;
	int ncopy;
	int copyfrom[NCH];

	// member functions
	corr(int c, char *v[]);
	int getdata(datastruct *newdata, unsigned int Nnewdata);	// get next set of data
	string report();											// returns report string
	string reportbyround();										// returns report string
	int hist2();												// new 2-fold correlation
	int save2(char *fn, char ch_a, char ch_b);					// save histogram for 2
	int hhist();												// high-order correlation
	int save3(char *fn, char ch_a, char ch_b, char ch_c);	// save histogram for 3
	int save3b(char *fn, char ch_a, char ch_b, char ch_c);	// save normalized histogram for 3
	int save4(char *fn, char ch_a, char ch_b, char ch_c, char ch_d);	// save histogram for 4
	int save4b(char *fn, char ch_a, char ch_b, char ch_c, char ch_d);	// save normalized histogram for 4
	void savetimestamp(char *fn);							// save (channel, timestamp) array

	// functions introduced in v18: piecewise meas values
	int getg2();
	int gethg();
	int getg2avg();
	int gethgavg();
	int saveg2avg(char *fn, char ch_a, char ch_b);							// save g2
	int saveg3avg(char *fn, char ch_a, char ch_b, char ch_c);				// save g3
	int saveg4avg(char *fn, char ch_a, char ch_b, char ch_c, char ch_d);	// save g4

};

// for sort
bool operator<(const datastruct& a, const datastruct& b)
{
	return ((a.timestamp != b.timestamp) ? (a.timestamp < b.timestamp) : (a.channel < b.channel));
}

corr::corr(int c, char *v[])
{
	int i,j,k;

	// check argument numbers
	if (c<6)   {
		cout << "Syntax: ht2hist <ht2 filename> <outfile prefix> <# of pts per round>";
		cout << " <time window in sec> <hist2 bin size in sec (ex: 1e-9)>";
		cout << " [-h 3 <hist3 bin size in sec> <channel#> <channel#> <channel#>]";
		cout << " [-c <channel#> [<channel#> ...]]";
		cout << " [-o <offset1 in sec> ... <offset4> [<offset5> ...] (ex: 1e-10)]";
		cout << " [-r <round# to stop after>] (ex: 10)]";
		cout << " [-noheader]\n";
		exit (0);
	}

	// process arguments
/*	want3=0; want4=0;
	switch (c) {
		case 9:
			want3 = 1;
			break;
		case 10:
			want3 = 1;
			want4 = 1;
			break;
	}*/
	if (c > 5) {

		// required parameters or defaults
		i=1;
		ftt = v[i++];
		fpre = v[i++];
		ndatareq = atoi(v[i++]);
		w = atof(v[i++]);
		binsize2 = (long long)(atof(v[i++])/TCONV);
		cout << "binsize2=" << binsize2 << "\n";
		// defaults
		highorder = 2;
		for (j=0; j<NCHTOT; j++)
			ttoffset[j]=0;
		wantoffset = 0;
		wantstopafter = 0;
		noheader = 0;
		for (j=0; j<NCH; j++)
			copyfrom[j] = 0;
		wantcopy =0;

		// optional parameters
		cout << c << " cmdline args given.\n";
		while (i < c) {
			cout << i << "\n";
			if (v[i]==string("-h")) {			// high order correlation parameters
				i++;
				highorder = atoi(v[i++]);
				hbinsize = (long long)(atof(v[i++])/TCONV);
				//cout << "hbinsize = " << hbinsize << "\n";
				for (j=0; j<highorder; j++)
					ch[j] = CH2I(atoi(v[i++]));
			} else if (v[i]==string("-c")) {	// channel copy parameters
				i++; ncopy=0;
				while (i<c) {
					if (v[i][0] == '-') break;
					copyfrom[ncopy] = atoi(v[i++]);
					cout << "copyfrom[" << ncopy << "]=" << copyfrom[ncopy] << "\n";
					ncopy++;
				}
				wantcopy = 1;
				cout << "-c processed\n";
			} else if (v[i]==string("-o")) {	// offset parameters
				i++;
				for (j=0; j<NCHTOT; j++) {
					ttoffset[j] = (long long)(atof(v[i++])/TCONV);
					cout << "ttoffset[" << j << "]=" << ttoffset[j] << "\n";
				}
				wantoffset = 1;
			} else if (v[i]==string("-r")) {	// round # to stop after
				i++;
				stopafterround = atoi(v[i++]);
				cout << "stopafterround=" << stopafterround << "\n";
				wantstopafter = 1;
			} else if (v[i]==string("-noheader")) {	// no header
				noheader = 1;
				cout << "noheader option given!\n";
				i++;
			} else {
				cout << "what is " << v[i] << " (" << i << "-th of " << c << " arguments)?\n";
				exit(0);
			}
		}
		cout << "Cmdline parsing done.\n";

/*		// higher order parameters
		if (c>6)
			if (((highorder=atoi(v[6])) > 2) && (highorder<=NCHTOT)) {
				hbinsize = (long long)(atof(v[7])/TCONV);
				for (i=0; i<highorder; i++)
					ch[i] = I2CHINV(atoi(v[i+8]));
			} else
				highorder = 2;
*/	} else
		cout << "Invalid command line parameters!\n";

	// timing constants, etc
	wbin=(w/TCONV);
	nbin2=(long)(wbin/binsize2)+1;
	if (highorder > 2)
		hnbin = (long)(wbin/hbinsize)+1;
	ncarryover = 0;
	for (i=0; i<NCHTOT; i++) ctot[i] = 0;
	ttot = 0;
	hnbintot = pow((double)(2*hnbin-3),(double)highorder-1);		// total number of bins in highorder histogram

	// tt data
	data = new datastruct [ndatareq*(1+wantcopy)];
	datacopy = new datastruct [ndatareq];

	// initialize 2nd order correlation variables
	for (i=0; i<NCHTOT; i++)
		for (j=0; j<NCHTOT; j++) {
			c2[i][j] = new unsigned int [2*nbin2-1];
			c2byround[i][j] = new unsigned int [2*nbin2-1];
			g2byround[i][j] = new double [2*nbin2-1];
			g2avg[i][j] = new double [2*nbin2-1];
			for (k=0; k<2*nbin2-1; k++) {
				c2[i][j][k]=0;
				c2byround[i][j][k]=0;
				g2byround[i][j][k]=0;
				g2avg[i][j][k]=0;
			}
		}
	cout << "Number of bins for 2-fold = " << 2*nbin2-3 << "\n";		// -N..0..+N and exclude boundaries

	// initialize higher order correlation variable
	if (highorder > 2) {
		cout << "Number of bins for " << highorder << "-fold = " << 2*hnbin-3 << "^" << highorder-1 << "\n";
		hc=new unsigned int [(unsigned int)pow((double)2*hnbin-1,highorder-1)];
		hcbyround=new unsigned int [(unsigned int)pow((double)2*hnbin-1,highorder-1)];
		hgbyround=new double [(unsigned int)pow((double)2*hnbin-1,highorder-1)];
		hgavg=new double [(unsigned int)pow((double)2*hnbin-1,highorder-1)];
		if (hc==NULL) exit (1);
		if (highorder == 3)
			for (i=0; i<2*hnbin-1; i++)
				for (j=0; j<2*hnbin-1; j++) {
					hc[(2*hnbin-1)*i+j]=0;
					hcbyround[(2*hnbin-1)*i+j]=0;
					hgbyround[(2*hnbin-1)*i+j]=0;
					hgavg[(2*hnbin-1)*i+j]=0;
				}
		else if (highorder == 4)
			for (k=0; k<2*hnbin-1; k++)
				for (j=0; j<2*hnbin-1; j++)
					for (i=0; i<2*hnbin-1; i++) {
						hc[i+(2*hnbin-1)*j+(2*hnbin-1)*(2*hnbin-1)*k]=0;
						hcbyround[i+(2*hnbin-1)*j+(2*hnbin-1)*(2*hnbin-1)*k]=0;
						hgbyround[i+(2*hnbin-1)*j+(2*hnbin-1)*(2*hnbin-1)*k]=0;
						hgavg[i+(2*hnbin-1)*j+(2*hnbin-1)*(2*hnbin-1)*k]=0;
					}
		// allocate memory for high order histogramming
		ndet = new unsigned int [highorder-1];
		tb = new long*[highorder-1];
		for (i=0; i<highorder-1; i++)
			tb[i] = new long[ndatareq];
	}

	// get max offsets
	long long pttoffset[NCHTOT];
	for (i=0; i<NCHTOT; i++)
		pttoffset[i] = ttoffset[i];
	sort(pttoffset, pttoffset+NCHTOT);
	ttoffsetmin = pttoffset[0];
	ttoffsetmax = pttoffset[NCHTOT-1];
	cout << "offsets: min=" << ttoffsetmin << " max=" << ttoffsetmax << "\n";

	cout << "Meas class created.\n";

	// end of memory allocation

}

// get next set of data.
// Let's have this function process all the information but correlation
int corr::getdata(datastruct *newdata, unsigned int Nnewdata)
{
	int i,j,k;
	unsigned long long ndatacopy;

	// duplicate channels if copy requested
	if (wantcopy) {
		ndatacopy = 0;
		for (i=0; i<Nnewdata; i++)
			for (j=0; j<ncopy; j++)
				if (newdata[i].channel == copyfrom[j]) {
					datacopy[ndatacopy].channel = newdata[i].channel + NCH;		// assign new channel number
					datacopy[ndatacopy].timestamp = newdata[i].timestamp;		// same timestamp
					ndatacopy++;
				}
		cout << ndatacopy << "\n";
	} else
		ndatacopy = 0;

	// gather data
	memcpy(&data[0], &data[ndata-ncarryover], sizeof(*data)*ncarryover);		// move carryover data to the front
	memcpy(&data[ncarryover], &newdata[0], sizeof(*data)*Nnewdata);				// add new data
	if (wantcopy)
		memcpy(&data[ncarryover+Nnewdata], &datacopy[0], sizeof(*data)*ndatacopy);	// add data copy

	unsigned int prevndata = ndata;
	ndata = ncarryover + Nnewdata + ndatacopy;								// set total # of data
	//savetimestamp("ts.dat");
	cout << ndata << "\n";

	// add offsets and sort the tt data
	if (wantoffset || wantcopy) {
		for (i=ncarryover; i<ndata; i++)
			data[i].timestamp += ttoffset[data[i].channel];
		cout << "Sorting...\n";
		sort(data, data+ndata);
	}

	// set istart
	if (irep == 0) {		// first iteration
		istart=1;
//		while ((data[istart].timestamp - data[0].timestamp) < wbin)
		while ((data[istart].timestamp - data[0].timestamp) < (ttoffsetmax - ttoffsetmin + wbin))
			istart++;
		tstart0 = data[istart].timestamp;	// very start time (corrected by time window!)
	} else
		istart = ncarryover - (prevndata - iend) + 1;		// iend comes from the last set of data

	// set iend and ncarryover
	iend = ndata-2;
	ncarryover = 2;
	while ((data[ndata-1].timestamp - data[iend].timestamp) < (ttoffsetmax - ttoffsetmin + wbin)) {
		iend--;
		ncarryover++;
	}
	i = iend-1;
	while ((data[iend].timestamp - data[i--].timestamp) < wbin)
		ncarryover++;

/*	//debug
	istart=0;
	iend=ndata-1;
	ncarryover=0;
	if (irep==0)
		tstart0=data[istart].timestamp;*/
//	cout << "tstart0=" << tstart0 << " tstart=" << data[istart].timestamp << " t(end+1)=" << data[iend+1].timestamp
//		<< " istart=" << istart << " iend=" << iend <<  " ncarryover=" << ncarryover << "\n";

	// accumulated measurement time
	ttot = data[iend].timestamp - tstart0;

	// accumulated number of clicks on each channel
	for (i=istart; i<=iend; i++)
		ctot[CH2I(data[i].channel)]++;
	
	// normalizations for g2 or g3
	for (j=0; j<NCHTOT; j++)
		if (ctot[j]>0)
			for (k=j; k<NCHTOT; k++) 
				if (ctot[k]>0)
					norm2[j][k] = (double)ctot[j]*ctot[k]*binsize2/ttot;
	//if (highorder == 3)
	//	norm3 = (double)ctot[ch[0]]*ctot[ch[1]]*ctot[ch[2]]*hbinsize*hbinsize/ttot/ttot;
	if (highorder > 2) {
		hnorm = (double)ctot[ch[0]];
		for (i=1; i<highorder; i++)
			hnorm *= (double)ctot[ch[i]]*hbinsize/ttot;
	}

	//// added on 02/09/09.
	// measurement time of this round
	tbyround = data[iend].timestamp - data[istart].timestamp;

	// counts of this round
	for (i=0; i<NCHTOT; i++) cbyround[i] = 0;
	for (i=istart; i<=iend; i++)
		cbyround[CH2I(data[i].channel)]++;
	//// end of the addition

	//// addition in v18 ////
	// reset histogram-by-round
	for (i=0; i<NCHTOT; i++)
		for (j=0; j<NCHTOT; j++)
			for (k=0; k<2*nbin2-1; k++)
				c2byround[i][j][k]=0;
	if (highorder > 2)
		if (highorder == 3)
			for (i=0; i<2*hnbin-1; i++)
				for (j=0; j<2*hnbin-1; j++)
					hcbyround[(2*hnbin-1)*i+j]=0;
		else if (highorder == 4)
			for (k=0; k<2*hnbin-1; k++)
				for (j=0; j<2*hnbin-1; j++)
					for (i=0; i<2*hnbin-1; i++)
						hcbyround[i+(2*hnbin-1)*j+(2*hnbin-1)*(2*hnbin-1)*k]=0;

	// normalizations by round
	for (j=0; j<NCHTOT; j++)
		if (ctot[j]>0)
			for (k=j; k<NCHTOT; k++) 
				if (ctot[k]>0)
					norm2byround[j][k] = (double)cbyround[j]*cbyround[k]*binsize2/tbyround;
	if (highorder > 2) {
		hnormbyround = (double)cbyround[ch[0]];
		for (i=1; i<highorder; i++)
			hnormbyround *= (double)cbyround[ch[i]]*hbinsize/tbyround;
	}
	//// end addition in v18 ////

	cout << "ndata=" << ndata << " Nnewdata=" << Nnewdata << " istart=" << istart << " iend=" << iend << "\n";
	cout << "truncated time: front=" << (data[istart].timestamp-data[0].timestamp)*TCONV
		<< " end=" << (data[ndata-1].timestamp-data[iend].timestamp)*TCONV << "\n";
	cout << "total measurement time = " << ttot*TCONV << "\n";
	cout << "measurement time of this round including truncations = " << (data[ndata-1].timestamp-data[0].timestamp)*TCONV << "\n";
	cout << "rate in this round = " << cbyround[0]/(TCONV*tbyround) << " " << cbyround[1]/(TCONV*tbyround) << " "
		<< cbyround[2]/(TCONV*tbyround) << " " << cbyround[3]/(TCONV*tbyround) << " Hz\n";
	cout << "timestamp[istart-1]=" << data[istart-1].timestamp << "\ttimestamp[istart]=" << data[istart].timestamp
		<< "\ttimestamp[iend]=" << data[iend].timestamp << "\ttimestamp[iend+1]=" << data[iend+1].timestamp << "\n";
	return 0;
}

// histogram for 2.
// take care of every channel combination at once.
int corr::hist2()
{
	unsigned int i,j;
	unsigned long long dt, bin;
	double nrbin;
	int ch1, ch2;
	clock_t startclock,finishclock;
	
	startclock = clock();
	// reset c2totbyround
	for (i=0; i<NCHTOT; i++)
		for (j=0; j<NCHTOT; j++)
			c2totbyround[i][j] = 0;

	for (i=istart; i<=iend; i++) {
		ch1 = CH2I(data[i].channel);

		// scan positive time
		if (i<ndata) {
			j=i+1;
			while ((dt = data[j].timestamp - data[i].timestamp) <= wbin) {
				if (data[j].channel >= data[i].channel) {					// only interested in 00,01,02,03,11,12,13,22,23
					bin = nbin2 - 1 + (long)((double)dt/binsize2+0.5);
					nrbin = nbin2 - 1 + (double)dt/binsize2+0.5;			// non-rounded
					//if ((nrbin-(double)bin)>0) {							// discard a count at bin-boundaries
						ch2 = CH2I(data[j].channel);
						c2[ch1][ch2][bin]++;
						c2byround[ch1][ch2][bin]++;
						if ((bin>0) && (bin<2*nbin2-2))
							c2totbyround[ch1][ch2]++;
					//}
				}
				j++;
			}
		}

		// scan negative time
		if (i>0) {
			j=i-1;
			while ((dt = data[i].timestamp - data[j].timestamp) <= wbin) {
				if (data[j].channel >= data[i].channel) {					// only interested in 00,01,02,03,11,12,13,22,23
					bin = (long)(nbin2 - 1 - (double)dt/binsize2 + 0.5);
					nrbin = nbin2 - 1 - (double)dt/binsize2+0.5;			// non-rounded
					//if ((nrbin-(double)bin)>0) {							// discard a count at bin-boundaries
						ch2 = CH2I(data[j].channel);
						c2[ch1][ch2][bin]++;
						c2byround[ch1][ch2][bin]++;
						if ((bin>0) && (bin<2*nbin2-2))
							c2totbyround[ch1][ch2]++;
					//}
				}
				j--;
			}
		}

		// show the progress every 10 sec.
		finishclock = clock();
		//cout << "| " << startclock << " " << finishclock << " | ";
		if (((double(finishclock)-double(startclock))/CLOCKS_PER_SEC) > 5) {
			cout << i << "," << (int)c2[0][1][0] << "  ";
			startclock = clock();
		}
	}
	return 0;
}

// high order correlation
int corr::hhist()
{
	int i,j,k,l;
	int dt, bin;
	double nrbin;
/*	long **tb;
	unsigned int *ndet;
*/	clock_t startclock,finishclock;
	
	startclock = clock();
	hctotbyround=0;		// count the total

/*	ndet = new unsigned int [highorder-1];
	tb = new long*[highorder-1];
	for (i=0; i<highorder-1; i++)
		tb[i] = new long[ndata];
*/
	for (i=istart; i<=iend; i++) {
		if (data[i].channel==ch[0]) {
			// reset # of detections
			for (j=0; j<highorder-1; j++)
				ndet[j] = 0;
			
			// scan for negative
			j=i-1;
			//j=i;
			while (((dt=data[i].timestamp-data[j].timestamp) <= wbin) && (j>0)) {
				k=1;
				do
				if (data[j].channel==ch[k]) {
					//cout << "other channel detected!\n";
					bin = (long)(hnbin - 1 - (double)dt/hbinsize+0.5);
					nrbin = hnbin - 1 - (double)dt/hbinsize+0.5;
					if ((bin>0) && (bin<(2*hnbin-2))) {		// exclude tw-boundary counts
						//if ((nrbin - (double)bin) > 0)		// exclude bin-boundary counts
							*(tb[k-1]+ndet[k-1]++) = bin-1;
							//k = highorder;	// cheat and get out of the loop
					}
				}
				while (++k < highorder);
				j--;
			}

			// scan for positive
			j=i+1;
			//j=i;
			while (((dt=data[j].timestamp-data[i].timestamp) <= wbin) && (j<ndata)) {
				k=1;
				do
				if (data[j].channel==ch[k]) {
					bin = (long)(hnbin - 1 + (double)dt/hbinsize+0.5);
					nrbin = hnbin - 1 + (double)dt/hbinsize+0.5;
					if ((bin>0) && (bin<(2*hnbin-2))) { 		// exclude boundary counts
						//if ((nrbin - (double)bin) > 0)		// exclude bin-boundary counts
							*(tb[k-1]+ndet[k-1]++) = bin-1;
							//k = highorder;	// cheat and get out of the loop
					}
				}
				while (++k < highorder);
				j++;
			}

			// calculate histogram 3
			if (highorder==3)
				for (k=0; k<ndet[1]; k++)
					for (j=0; j<ndet[0]; j++) {
						hc[tb[0][j] + tb[1][k]*(2*hnbin-3)]++;
						hcbyround[tb[0][j] + tb[1][k]*(2*hnbin-3)]++;
						hctotbyround++;		// count the total
					}
			// histogram 4
			else if (highorder==4)
				for (l=0; l<ndet[2]; l++)
					for (k=0; k<ndet[1]; k++)
						for (j=0; j<ndet[0]; j++) {
							hc[tb[0][j] + tb[1][k]*(2*hnbin-3) + tb[2][l]*(2*hnbin-3)*(2*hnbin-3)]++;
							hcbyround[tb[0][j] + tb[1][k]*(2*hnbin-3) + tb[2][l]*(2*hnbin-3)*(2*hnbin-3)]++;
							hctotbyround++;		// count the total
						}
		}

		// show the progress every 10 sec.
		finishclock = clock();
		//cout << "| " << startclock << " " << finishclock << " | ";
		if (((double(finishclock)-double(startclock))/CLOCKS_PER_SEC) > 5) {
			//debug
			//for (i=0; i<hnbin; i++)
			//	cout << hc[i*hnbin] << " ";

			cout << i << "," << hc[0] << "  ";
			startclock = clock();
		}
	}
/*	for (i=0; i++; i<highorder)
		delete[] tb[i];
	delete[] tb;
	delete[] ndet;
*/
	return 0;
}

// save the summary
string corr::report()
{
	int i,j,k;
	stringstream ssin("",ios_base::out);

	ssin << "Report =======================\n"
		<< "ht2hist version: " << VERSION << "\n"
		<< "ht2 file: "<< ftt << "\n"
		<< "Rounds: " << irep+1 << "\n"
		<< "Net measurement time (s):  " << TCONV*ttot << "\n"
		<< "Time window (s): " << TCONV*wbin << "\n"
		<< "Bin size (2nd order): " << binsize2*TCONV << "\n"
		<< "Number of bins (2nd order) = " << 2*nbin2-3 << "\n";		// -N..0..+N and exclude boundaries
	if (highorder > 2) {
		ssin << "Bin size (" << highorder << "-th order) = " << hbinsize*TCONV << "\n"
			<< "Number of bins (" << highorder << "-th order) = " << 2*hnbin-3 << "^" << highorder << "\n";
	}
	
	ssin << "Copied channel: ";
	for (i=0; i<ncopy; i++)
		ssin << copyfrom[i] << "  ";
	ssin << "\n";

	for (i=0; i<NCHTOT; i++)
		ssin << "Time offset " << i << " (ps): " << ttoffset[i] << "\n";

	ssin << "Count_ch0: " << ctot[0] << "\n"
		<< "Count_ch1: " << ctot[1] << "\n"
		<< "Count_ch2: " << ctot[2] << "\n"
		<< "Count_ch3: " << ctot[3] << "\n"
		<< "Count_ch4: " << ctot[4] << "\n"
		<< "Count_total: " << ctot[0]+ctot[1]+ctot[2]+ctot[3]+ctot[4] << "\n"
		<< "Rate_ch0: " << ctot[0]/(TCONV*ttot) << "\n"
		<< "Rate_ch1: " << ctot[1]/(TCONV*ttot) << "\n"
		<< "Rate_ch2: " << ctot[2]/(TCONV*ttot) << "\n"
		<< "Rate_ch3: " << ctot[3]/(TCONV*ttot) << "\n"
		<< "Rate_ch4: " << ctot[4]/(TCONV*ttot) << "\n"
		<< "Rate_total: " << (ctot[0]+ctot[1]+ctot[2]+ctot[3]+ctot[4])/(TCONV*ttot) << "\n"
		<< "\n"
		<< "Normalization for 2:\n";
	for (j=0; j<NCHTOT; j++)
		if (ctot[j]>0)
			for (k=j; k<NCHTOT; k++) 
				if (ctot[k]>0)
					ssin << "   n" << I2CH(j) << I2CH(k) << " = " << norm2[j][k] << "\n";
	if (highorder > 2) {
		ssin << "Normalization for " << highorder << ":\n" << "   n";
		for (i=0; i<highorder; i++)
			ssin << I2CH(ch[i]);
		ssin << " = " << hnorm << "\n";
	}
	ssin << "End of report ================\n";
	return ssin.str();
}

// save a line report of every round
string corr::reportbyround()
{
	int i,j,k;
	unsigned long long c2tott;
	unsigned long binzdel;

	stringstream ssin("",ios_base::out);

	ssin << irep+1 << "\t"
		<< data[istart].timestamp*TCONV << "\t"
		<< data[iend].timestamp*TCONV << "\t"
		<< tbyround*TCONV << "\t"
		<< cbyround[0] << "\t"
		<< cbyround[1] << "\t"
		<< cbyround[2] << "\t"
		<< cbyround[3] << "\t"
		<< cbyround[4] << "\t";

	//hist2
	for (j=0; j<NCH-1; j++)
		for (k=j+1; k<NCH; k++) {
			ssin << c2[j][k][nbin2-1] <<"\t";			// accumulated count at zero delay
			ssin << (double)c2totbyround[j][k]/(2*nbin2-3) << "\t";			// average count in this round
		}

	//hist3 (or 4)
	if (highorder > 2) {
		binzdel = 0;
		for (i=0; i<(highorder-1); i++)
			binzdel += pow((double)2*hnbin-3,i)*(hnbin-2);
		ssin << hc[binzdel] << "\t";		// accumulated count at zero delay
		ssin << hctotbyround/pow((double)(2*hnbin-3),(highorder-1));		// average count in this round
	}
	ssin << "\n";
	return ssin.str();
}

// save 2nd order correlation c2
int corr::save2(char *fn, char ch_a, char ch_b)
{
	int i;
	ofstream f;

	f.open(fn);
	for (i=1; i<2*nbin2-2; i++) {			// i=0 and 2*nbin2-2 are boundary bins with half counts; throw them out!
		f << (i-nbin2+1)*binsize2*TCONV << "\t" << c2[ch_a][ch_b][i] << "\t"
			<< c2[ch_a][ch_b][i]/norm2[ch_a][ch_b] << "\t"			// normalized and accumulated
			<< g2avg[ch_a][ch_b][i] << "\n";				// g2avg
			//<< c2[ch_a][ch_b][i]/c0*tlim/(tlim-binsize2*abs(i-nbin2+1)) << "\n";		// normalized and weighted
	}
	f.close();
	cout << fn << " ";
	return 0;
}

// save 3rd order correlation c3
int corr::save3(char *fn, char ch_a, char ch_b, char ch_c)
{
	int i,j;
	ofstream f;
	f.open(fn);	
	
	// normalized matrix form
	for (i=0; i<2*hnbin-3; i++) {
		for (j=0; j<2*hnbin-3; j++)
			f << hc[(2*hnbin-3)*i+j] << "\t";
		f << "\n";
	}
	f.close();
	cout << fn << " ";
	return 0;
}

// save normalized hist3
int corr::save3b(char *fn, char ch_a, char ch_b, char ch_c)
{
	int i,j;
	ofstream f;
	f.open(fn);	
	
	for (i=0; i<2*hnbin-3; i++) {
		for (j=0; j<2*hnbin-3; j++)
			f << (double)(hc[(2*hnbin-3)*i+j])/hnorm << "\t";
		f << "\n";
	}
	f.close();
	cout << fn << "\t";
	return 0;
}

// save 4th order correlation c3
int corr::save4(char *fn, char ch_a, char ch_b, char ch_c, char ch_d)
{
	int i,j,k;
	ofstream f;

	f.open(fn, ios::out | ios::binary);
	unsigned int ndatasave = (2*hnbin-3)*(2*hnbin-3)*(2*hnbin-3);
	cout << "sizeof(*hc)=" << sizeof(*hc) << " ndatasave=" << ndatasave << "\n";
	f.write((char *)hc, sizeof(*hc)*ndatasave);
	
	/*	// normalized matrix form
	for (k=1; k<2*hnbin-2; k++) {			// exclude boundary bins
		for (j=1; j<2*hnbin-2; j++)
			for (i=1; i<2*hnbin-2; i++) 
				f << hc[i + (2*hnbin-1)*j + (2*hnbin-1)*(2*hnbin-1)*k] << "\n";
		//f << "\n";
	}*/

	f.close();
	cout << fn << " ";
	return 0;
}

// save normalized 4th order correlation c3
int corr::save4b(char *fn, char ch_a, char ch_b, char ch_c, char ch_d)
{
	int i,j,k;
	ofstream f;
	f.open(fn);	
	
	// normalized matrix form
	for (k=1; k<2*hnbin-2; k++) {			// exclude boundary bins
		for (j=1; j<2*hnbin-2; j++)
			for (i=1; i<2*hnbin-2; i++)
				f << hc[i + (2*hnbin-1)*j + (2*hnbin-1)*(2*hnbin-1)*k]/hnorm << "\n";
		//f << "\n";
	}
	f.close();
	cout << fn << " ";
	return 0;
}

void corr::savetimestamp(char *fn)
{
	int i;
	ofstream f;
	f.open(fn);
	for (i=0; i<ndata; i++) {
		f << (int)data[i].channel << "\t" << data[i].timestamp << "\t" << "\n";
	}
	f.close();
	cout << " " << fn << " ";
}

//////////////////////////
// new functions in v18
//////////////////////////
int corr::getg2()
{
	int ch_a,ch_b,i;

	for (ch_a=0; ch_a<4; ch_a++)					// channel A
		if (ctot[ch_a]>0)							// if ch A was on
			for (ch_b=ch_a; ch_b<4; ch_b++)			// channel B
				if (ctot[ch_b]>0)					// if ch B was on
					for (i=0; i<2*nbin2-1; i++) {
						//cout << "in a loop of g2byround.\n";
						g2byround[ch_a][ch_b][i] = c2byround[ch_a][ch_b][i]/norm2byround[ch_a][ch_b];	// g2 of A and B
					}
	return 0;
}

int corr::gethg()
{
	int i;

	for (i=0; i<hnbintot; i++)
		hgbyround[i] = (double)(hcbyround[i])/hnormbyround;
	return 0;
}

int corr::getg2avg()
{
	int ch_a,ch_b,i;

	for (ch_a=0; ch_a<NCH; ch_a++)					// channel A
		if (ctot[ch_a]>0)							// if ch A was on
			for (ch_b=ch_a; ch_b<NCH; ch_b++)			// channel B
				if (ctot[ch_b]>0)					// if ch B was on
					for (i=0; i<2*nbin2-1; i++)
						g2avg[ch_a][ch_b][i] = (g2avg[ch_a][ch_b][i]*irep + g2byround[ch_a][ch_b][i])/(irep+1);	// update g2avg
	return 0;
}

int corr::gethgavg()
{
	int i;

	for (i=0; i<hnbintot; i++)
		hgavg[i] = (hgavg[i]*irep + hgbyround[i])/(irep+1);
	return 0;
}

int corr::saveg2avg(char *fn, char ch_a, char ch_b)							// save g2
{
	int i;
	ofstream f;

	f.open(fn);
	for (i=1; i<2*nbin2-2; i++) {			// i=0 and 2*nbin2-2 are boundary bins with half counts; throw them out!
		f << (i-nbin2+1)*binsize2*TCONV << "\t" << g2avg[ch_a][ch_b][i] << "\n";		// normalized
			//<< c2[ch_a][ch_b][i]/c0*tlim/(tlim-binsize2*abs(i-nbin2+1)) << "\n";		// normalized and weighted
	}
	f.close();
	cout << fn << " ";
	return 0;
}

int corr::saveg3avg(char *fn, char ch_a, char ch_b, char ch_c)				// save g3
{
	int i,j;
	ofstream f;
	f.open(fn);	
	
	for (i=0; i<2*hnbin-3; i++) {
		for (j=0; j<2*hnbin-3; j++)
			f << hgavg[(2*hnbin-3)*i+j] << "\t";
		f << "\n";
	}
	f.close();
	cout << fn << "\t";
	return 0;
}

int corr::saveg4avg(char *fn, char ch_a, char ch_b, char ch_c, char ch_d)	// save g4
{
	int i,j,k;
	ofstream f;

	f.open(fn, ios::out | ios::binary);
	unsigned int ndatasave = (2*hnbin-3)*(2*hnbin-3)*(2*hnbin-3);
	cout << "sizeof(*hc)=" << sizeof(*hgavg) << " ndatasave=" << ndatasave << "\n";
	f.write((char *)hgavg, sizeof(*hgavg)*ndatasave);
	f.close();
	cout << fn << " ";
	return 0;
}


////////////////////////////
// Main
////////////////////////////
int main(int argc, char *argv[])
{
    int i, j, k;
	unsigned int nreq;
	char fn[100];
	clock_t startclock,finishclock;

	cout << "Program start\n";
	cout << "Version " << VERSION << "\n";

	// main class
	corr meas(argc, argv);

	// handle data file (ht2)
	ht2 ttdatasrc;
	// open and get data file information
	if (ttdatasrc.openfile(argv[1]) != 0) {
		printf("error opening file!\n");
		return -1;
	}
	if (meas.noheader ==1)
		cout << "No header reading.\n\n";
	else {
		cout << "Header reading started.\n";
		ttdatasrc.ndatareq = meas.ndatareq;
		ttdatasrc.getHdr();
		cout << "Header reading done.\n\n";
	}

	startclock = clock();
	meas.irep = 0;
	while (ttdatasrc.chkeof() == 0) {

		// progress indicator
		if ((meas.irep%1)==0)
			cout << "Round " << meas.irep+1 << "\n";
		if (_kbhit())
			if (_getch()=='q')
				break;
		if (meas.wantstopafter==1)
			if (meas.irep+1 > meas.stopafterround)
				break;

		// fetch tt data
		nreq = meas.ndatareq - meas.ncarryover;
		ttdatasrc.getTTData(nreq);
		meas.getdata((datastruct *)ttdatasrc.ttdata, ttdatasrc.ttdatasize);	// transfer data to the main object
		
		// histogram 2nd order correlation
		cout << "Hist-2 started.\n";
		meas.hist2();
		cout << "Hist-2 done.\n";
		cout << "count2=" << meas.c2[0][0][0] << "\n";
		meas.getg2();
		meas.getg2avg();

		// histogram higher order correlation
		if (meas.highorder > 2) {
			cout << "Hist-" << meas.highorder << " of ";
			for (i=0; i<meas.highorder; i++)
				cout << I2CH(meas.ch[i]) << " ";
			cout << "started.\n";
			meas.hhist();
			cout << "Hist-" << meas.highorder << " done.\n";
			cout << "count" << meas.highorder << "=" << meas.hc[0] << "\n";
			meas.gethg();
			meas.gethgavg();
		}

		//cout << "hnorm = " << meas.hnorm << "\n";
		//cout << "norm2byround = " << meas.norm2byround[1][2] << "\n";
		//cout << "hnormbyround = " << meas.hnormbyround << "\n";
		cout << "g2 (g3, g4) calculation done.\n";

		// save the report
		ofstream f;
		sprintf(fn,"%sreport.txt",meas.fpre);	// compose filename
		f.open(fn);
		f << meas.report();
		f.close();
		// save histogram for 2
		for (j=0; j<NCHTOT; j++)
			if (meas.ctot[j]>0)
				for (k=j; k<NCHTOT; k++) 
					if (meas.ctot[k]>0) {
						//hist2(j,k);
						sprintf(fn,"%s%d%d.dat",meas.fpre,I2CH(j),I2CH(k));	// compose filename
						meas.save2(fn,j,k);
						//sprintf(fn,"%s%d%davg.dat",meas.fpre,I2CH(j),I2CH(k));	// compose filename
						//meas.saveg2avg(fn,j,k);
					}
		cout << "\n";

		// save timestamp data for debugging
		//sprintf(fn,"%s%TS.dat",meas.fpre);
		//meas.savetimestamp(fn);

		// save histogram for 3
		if (meas.highorder == 3) {
			sprintf(fn,"%s%d%d%d.dat",meas.fpre,I2CH(meas.ch[0]),I2CH(meas.ch[1]),I2CH(meas.ch[2]));
			meas.save3(fn,meas.ch[0],meas.ch[1],meas.ch[2]);
			sprintf(fn,"%s%d%d%dnorm.dat",meas.fpre,I2CH(meas.ch[0]),I2CH(meas.ch[1]),I2CH(meas.ch[2]));
			meas.save3b(fn,meas.ch[0],meas.ch[1],meas.ch[2]);
			sprintf(fn,"%s%d%d%davg.dat",meas.fpre,I2CH(meas.ch[0]),I2CH(meas.ch[1]),I2CH(meas.ch[2]));
			meas.saveg3avg(fn,meas.ch[0],meas.ch[1],meas.ch[2]);
			cout << "\n";
		}

		// save histogram for 4
		if (meas.highorder == 4) {
			sprintf(fn,"%s%d%d%d%d.dat",meas.fpre,I2CH(meas.ch[0]),I2CH(meas.ch[1]),I2CH(meas.ch[2]),I2CH(meas.ch[3]));
			meas.save4(fn,meas.ch[0],meas.ch[1],meas.ch[2],meas.ch[3]);
			//sprintf(fn,"%s%d%d%d%dnorm.dat",meas.fpre,I2CH(meas.ch[0]),I2CH(meas.ch[1]),I2CH(meas.ch[2]),I2CH(meas.ch[3]));
			//meas.save4b(fn,meas.ch[0],meas.ch[1],meas.ch[2],meas.ch[3]);
			sprintf(fn,"%s%d%d%d%davg.dat",meas.fpre,I2CH(meas.ch[0]),I2CH(meas.ch[1]),I2CH(meas.ch[2]),I2CH(meas.ch[3]));
			meas.saveg4avg(fn,meas.ch[0],meas.ch[1],meas.ch[2],meas.ch[3]);
			cout << "\n";
		}
		cout << "\n";

		// save information of this data set
		sprintf(fn,"%sbyround.txt",meas.fpre);	// compose filename
		if (meas.irep==0)
			f.open(fn);
		else
			f.open(fn,ios::app);
		f << meas.reportbyround();
		f.close();

		meas.irep++;
	}
	meas.irep--;

	// show report
	cout << meas.report();
	cout << "Done. Bye!\n";

	return 0;
}
