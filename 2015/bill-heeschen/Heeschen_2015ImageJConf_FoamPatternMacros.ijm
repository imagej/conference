/*  Foam Pattern macros
Written by Bill Heeschen
The Dow Chemical Company
This code is for demonstration of the principals for analyzing foam patterns.
No warranty of applicability to any images is made, express or implied, by the author
or The Dow Chemical Company.

For presentation to the 2015 ImageJ Conference, Madison, WI

Assumes all necessary preprocessing of the image has been accomplished:
crop/rotate/center pattern/background correction, etc.

Requires the following ImageJ plugins:
	** Find_Ridges
		[Find_Ridges is from Bob Dougherty, OptiNav, Inc. 4/30/2014.
		Corrected coefficient of mixed partial derivative formula 5/3/2014]
	** Morphology  package from Gabriel Landini

*/
var tab="\t", sp=" ", nl="\n", uScore="_", mt="Empty";//abbreviations
var trueStr="true", falseStr="false";
var twoPI=2*PI;
var initDone=false;

var outlierRadius=12;//For outlier rejection routine
var outlierThresh=6;

var useVarRidges=true;
var varianceRad=16;
var ridgeRadVar=4;//for doRidges()
var ridgeThinVar=20;

var useLightRidges=true;
var useDarkRidges=true;
var ridgeRad=4;
var ridgeThin=20;

var meanRadius=10;//for flatten() routine
var ridgePix_A=75;//# of pixels in any skeleton to keep as a ridge candidate.
var ridgePix_B=40;//skeleton branches smaller than this will be pruned
var keeperPix=2*ridgePix_A;//If a branch is this long, keep it, even if it is part of a "node" structure

var refX = 0, refY = 0; refDist=0;//these are the reference points for the line arrays. 
var normArea=1;
var ANGLE_LUT=newArray(1);//placeholder for actual array
var anglePix=7;//number of pixels to use in the pixel-by-pixel angle calculation

var refAreaPix=1000000;//divide image size by this to get the relative area of the image

var variancePre="Var_";
var mainPre="Main_";
var invertPre="Invert_";
var ridgePre="Ridge_";
var ROIvariancePre="ROI_Variance_";
var ROIbrightPre="ROI_Bright_";
var ROIdarkPre="ROI_Dark_";
var colorVariance="green";
var colorBright="cyan";
var colorDark="blue";

//for profiles
var firstEndIndex=-1;
var secondEndIndex=-1;
var selectionXarr=newArray(1);
var selectionYarr=newArray(1);

//ROI characterization values
var midOffsetThresh = 0.08;//
var dA_ROIend_line=180/16;//PI/16   ROI_LocAngle and endsAngle must differ by less than this for the ROI to be a radial line
var radiusFactor=2.0;//allows angle range to "roll" toward long axis of image as distance from center grows.
var minARF_line=4.0;//Feret-based aspect ratio must be greater than this to be a true line 
var maxCurl_line=1.25;//curl calculated from end-to-end vs. length
var lineType_ScoreMax=3;//out of three - this is a rough guess
var lineType_ScoreThresh=2.8;//out of three - this is a rough guess
var minDevLengthRatio = 0.05;  //move this to globals
var devLengthRatioThresh=0.12;//  has to look like a chevron!
var dA_DevAngleDiff=180/16;//PI/16   ROI_LocAngle and devAngle must differ by less than this for the ROI to be a radial line
var chevron_ScoreMax=1;
var chevron_ScoreThresh=0.98;

//Reporting
var currDir="";
var inset=5, logW=380, logH=250;//default size of the Log window
var reportHeader="File"+tab+"LineCount"+tab+"ChevronCount"+tab+"SummedLength"+tab+"RidgeIndex";
var liveTableName="Measurements";
var liveTableNameB="["+liveTableName+"]";
var liveTableWide=600;
var liveTableHigh=300;

//for recording preferences
var STDPREF="foam_pat_pref.";
var keyA="typeA";//keys for IJ.prefs   STDPREF+key 
var typeAarr=newArray(1);//only need this as a global if actually USING the array as a global
var TYPE_A_STR= buildParamString(sp,false);//build this as the default for "Factory Settings" restore before the globals get replaced from IJ_Prefs.txt

macro "AutoRun" {
	run("Misc...", "divide=Infinity use run enhanced");//make sure the arrow cursor gets chosen.  Many of these images are very difficult to see the cross cursor!
}//end macro AutoRun

macro "Set Up Analysis [F1]"{
	if(!initDone) doInit();
	else doSetupDialog();//doInit already does the doSetupDialog
}//end macro  Set Up Analysis[] 

function doInit(){
	if(anglePix%2!=1)anglePix++;//if not odd, make it odd by making it larger
	ANGLE_LUT=makeAngleLUT(anglePix);
	prefReadWrite(true);//loads prefs from %USER\.imagej\IJ_Prefs.txt
	doSetupDialog();
	setUpReport();
	initDone=true;
}//end doInit()

macro "Ridges [F2]"{
	if(!initDone) doInit();
	if(nImages()<1) open();
	currID=getImageID();
	currDir=getDirectory("image");
	currName=getTitle();
	currUnit=""; currPixWide=0; currPixHigh=0;//keep these around in case I need to report in real-world units.
	getPixelSize(currUnit, currPixWide, currPixHigh);
	mainName=mainPre+currName;
	invertName=invertPre+currName;
	varianceName=variancePre+currName;
	run("8-bit");//just in case
	run("Scale...", "x=0.5 y=0.5 width=1413 height=702 interpolation=Bilinear average create title="+mainName);
	mainID=getImageID();
	run("Set Scale...", "distance=0 known=0 pixel=1 unit=pixel");
	mainW=getWidth();
	mainH=getHeight();
	normArea=mainW*mainH/refAreaPix;
	refX=round(mainW/2);
	refY=round(mainH/2);
	refDist=refX;	if(refY>refX)refDist=refY;
	run("Duplicate...", "title=working");
	workID=getImageID();
	run("Enhance Contrast", "saturated=0.1");
	run("Apply LUT");
	run("Mean...", "radius=5");
	flatten(workID,outlierRadius,outlierThresh, meanRadius);//ridges/valleys now pre-conditioned and contained in workID
	roiManager("reset");
	showStatus("Finding Ridges");
	if(useVarRidges){
		selectImage(workID);
		run("Duplicate...", "title="+varianceName);//get all ridges based on variance boundaries
		run("32-bit");
		varianceID=getImageID();
		run("Variance...", "radius="+d2s(varianceRad,0));
		run("8-bit");//variance image is now ready for ridge-finding
		returnImageID=doRidges(varianceID, ridgeRadVar, ridgeThinVar, false);
		roiManager("Deselect");
		roiManager("Save", currDir+ROIvariancePre+currName+".zip");
		roiManager("Set Color", colorVariance);
		roiManager("Set Line Width", 0);
		selectImage(mainID);
		roiManager("Show all without labels");
		run("From ROI Manager");//add variance ridges to Overlay
		roiManager("reset");
		closeByID(varianceID);
	}
	if(useDarkRidges){
		selectImage(workID);
		run("Duplicate...", "title="+invertName);
		run("Invert");
		invertID=getImageID();
		returnImageID=doRidges(invertID, ridgeRad, ridgeThin, false);
		roiManager("Deselect");
		roiManager("Save", currDir+ROIdarkPre+currName+".zip");
		roiManager("Set Color", colorDark);
		roiManager("Set Line Width", 0);
		selectImage(mainID);
		roiManager("Show all without labels");
		run("From ROI Manager");//add dark ridges to Overlay
		roiManager("reset");
		closeByID(invertID);
	}
	if(useLightRidges){
		returnImageID=doRidges(workID, ridgeRad, ridgeThin, false);
		roiManager("Deselect");
		roiManager("Save", currDir+ROIbrightPre+currName+".zip");
		roiManager("Set Color", colorBright);
		roiManager("Set Line Width", 0);
		selectImage(mainID);
		roiManager("Show all without labels");
		run("From ROI Manager");//add bright ridges to Overlay
		roiManager("reset");
	}
	countsForRidgeIndex=0;//if using light and dark ridges, then calculate mean/Std.Dev. of corresponding grays to report Ridge Index value
	if(useLightRidges) {
		roiManager("Open", currDir+ROIbrightPre+currName+".zip");
	}
	if(useDarkRidges) {
		roiManager("Open", currDir+ROIdarkPre+currName+".zip");
		if(useLightRidges) countsForRidgeIndex=roiManager("count");
	}
	if(useVarRidges){
		roiManager("Open", currDir+ROIvariancePre+currName+".zip");
	}
	IJ.deleteRows(0, nResults-1);//clear Results table
	run("Set Measurements...", "area mean standard modal min centroid perimeter fit feret's median redirect="+mainName+" decimal=3");
	roiManager("Deselect");
	roiManager("Measure");//with no ROIs selected, measures all ROIs
	numROIs=roiManager("count");
	showStatus("evaluating ridges");
	checkRefXY=false;
	typeArr=newArray(numROIs);
	lineCount=0;
	chevronCount=0;
	summedLength=0;
	meanArray=newArray(countsForRidgeIndex);
	ridgeIndex=0;
	ridgeFactor=NaN;
	if((!isNaN(refX))&&(!isNaN(refY))&&(refX>=0)&&(refY>=0)) checkRefXY=true;//if refX and refY >=0, then adjust X,Y arrays so that point closest to refX, refY is start of array
	for(ii=0;ii<numROIs;ii++){
		outFlag=findEnds(ii);
		if(!outFlag) {//could not find ends
			roiType=-1;
		}
		else {//ends found, so measure ROI for angles, etc.
			selectionXarr=Array.slice(selectionXarr,firstEndIndex,secondEndIndex+1);
			selectionYarr=Array.slice(selectionYarr,firstEndIndex,secondEndIndex+1);
			if(checkRefXY) {//if refX and refY >=0, then adjust X,Y arrays so that point closest to refX, refY is start of array
				lastInd=selectionXarr.length-1;
				dXone=selectionXarr[0]-refX;
				dYone=selectionYarr[0]-refY;
				dXtwo=selectionXarr[lastInd]-refX;
				dYtwo=selectionYarr[lastInd]-refY;
				if((dXone*dXone+dYone*dYone)>(dXtwo*dXtwo+dYtwo*dYtwo)){//assume first index is closer.  Only change if second index is closer.
					selectionXarr=Array.reverse(selectionXarr);
					selectionYarr=Array.reverse(selectionYarr);
				}
				//now correct ellipse fit angle to match reference coordinates
				//Note that all ellipse calculations result in directions of 0<=Angle<180
				currX=getResult("X",ii);
				currY=getResult("Y",ii);
				currAngle=getResult("Angle",ii);
				if(currX>=refX) { if(currAngle>90) setResult("Angle",ii,currAngle-180); }
				else if(currX<refX){ if(currAngle<90) setResult("Angle",ii,currAngle-180); }
			}
			tAngles=getAngleArr(selectionXarr,selectionYarr,anglePix);//anglePix is the number of pixels to use around the center pixel
			outFlag=getROIvalues(ii,tAngles,refDist);//this routine relies on the Results table generated from the roiManager("Measure") call
			roiType=assessROI(ii);//relies on Results table being populated by getROIvalues
		}
		setResult("ROItype",ii,roiType);
		if(roiType>=0){
			summedLength+=getResult("ROIlength",ii);
			if(roiType==1) lineCount++;
			else if(roiType==0) chevronCount++;
			if(ii<countsForRidgeIndex){//if true: this feature is from light or dark ridge, so include it in the RidgeIndex calculation
				meanArray[ridgeIndex]=getResult("Mean",ii);
				ridgeIndex++;//this is the current count in the array
			}
		}
	}
	closeByID(workID);//now close this one
	for(ii=numROIs-1;ii>=0;ii--){
		if(getResult("ROItype",ii)<0){
			IJ.deleteRows(ii,ii);
			roiManager("select",ii);
			roiManager("Delete");
		}
	}
	numROIs=roiManager("count");
	if(countsForRidgeIndex>0){
		meanArray=Array.trim(meanArray,ridgeIndex);//ridgeIndex is the current count
		Array.getStatistics(meanArray,ridgeMin, ridgeMax, ridgeMean, ridgeStDev);
		ridgeFactor=ridgeStDev*summedLength;
	}
	reportStr=currName+tab+d2s(lineCount/normArea,3)+tab+d2s(chevronCount/normArea,3)+tab+d2s(summedLength/normArea,3);
	if(!isNaN(ridgeFactor))reportStr=reportStr+tab+d2s(ridgeFactor/normArea,1);
	if(!isOpen(liveTableName)) setUpReport();
	print(liveTableNameB,reportStr);
	run("Hide Overlay");
	roiManager("Show All");
	showStatus("Ridges done");
}//end macro Ridges[]

macro "Reset macro parameters to factory default [F8]"{
	if(getBoolean("This will restore parameters"+nl+"to hard-coded defaults and update ImageJ Preferences"+nl+"to these values.  Continue?")){
		call("ij.Prefs.set",STDPREF+keyA, TYPE_A_STR);
		call("ij.Prefs.savePreferences");//do this explicitly for safety in case of crash
		result=extractParamString(TYPE_A_STR,sp);//return true if successful, false if it failed
		if(result) showMessage("Parameters restored to factory default.");
		else {
			showMessage("Problem restoring parameters - use dialog");
			doSetupDialog();
		}
	}
	else showMessage("No reset. Current parameters retained.");
}//end macro Reset macro parameters...[]

function doSetupDialog(){
	//add a configuration dialog here which includes an option to update the reference profile
	Dialog.create("Foam Pattern Analysis Setup");
	Dialog.addMessage("*** Ridge-finding parameters ***");
	Dialog.addNumber("Grayscale Outlier Radius",outlierRadius);//radius for outlier rejection routine
	Dialog.addNumber("Outlier Threshold", outlierThresh);//gray rejection range
	Dialog.addNumber("Smoothing Radius", meanRadius);//for flatten() routine
	Dialog.addCheckbox("Include Light ridges?", useLightRidges);
	Dialog.addCheckbox("Include Dark ridges?", useDarkRidges);
	Dialog.addCheckbox("Include Variance-based ridges?", useVarRidges);
	Dialog.addMessage("*** Parameters for light/dark ridges ***");
	Dialog.addNumber("Radius for doRidges routine", ridgeRad);
	Dialog.addNumber("Ridge thinning", ridgeThin);
	Dialog.addMessage("*** Parameters for variance-based ridges ***");
	Dialog.addNumber("Variance Radius", varianceRad);
	Dialog.addNumber("Radius for doRidges (variance)", ridgeRadVar);
	Dialog.addNumber("Ridge thinning (variance)", ridgeThinVar);
	Dialog.addMessage("*** Skeleton Parameters ***");
	Dialog.addNumber("Minimum pixels for a full skeleton", ridgePix_A);//# of pixels in a skeleton piece to keep as a ridge.
	Dialog.addNumber("Minimum pixels for a skeleton branches", ridgePix_B);//similar, but these are branches from the separated skeleton, so I'll be less stringent
	Dialog.addMessage("*** Feature classification parameters ***");
	Dialog.addNumber("Line - max curl", maxCurl_line);
	Dialog.addNumber("Line - min aspect", minARF_line);
	Dialog.addNumber("Line - max angle dev", dA_ROIend_line);
	Dialog.addNumber("Line - roll-off factor", radiusFactor);
	Dialog.addNumber("Chevron - min bulge ratio", devLengthRatioThresh);
	Dialog.addNumber("Chevron - max bulge angle dev", dA_DevAngleDiff);
	Dialog.addCheckbox("Print parameters to Log Window", false);
	Dialog.show();
//Read parameters
	//preprocess for ridges
	tReal=Dialog.getNumber();
	if(tReal>0) outlierRadius=tReal;//OK to be rational
	tInt=round(Dialog.getNumber());//gray rejection range
	if(tInt>0) outlierThresh=tInt;
	tReal=Dialog.getNumber();//for flatten() routine
	if(tReal>0) meanRadius=tReal;
	//find ridges
	useLightRidges=Dialog.getCheckbox();
	useDarkRidges=Dialog.getCheckbox();
	useVarRidges=Dialog.getCheckbox();
	tInt=round(Dialog.getNumber());
	if(tInt>0) ridgeRad=tInt;
	tInt=round(Dialog.getNumber());
	if(tInt>0) ridgeThin=tInt;
	tReal=Dialog.getNumber();
	if(tReal>0) varianceRad=tReal;//can be rational
	tInt=round(Dialog.getNumber());
	if(tInt>0) ridgeRadVar=tInt;
	tInt=round(Dialog.getNumber());
	if(tInt>0) ridgeThinVar=tInt;
	//process ridges
	tInt=round(Dialog.getNumber());//# of pixels in a skeleton piece to keep as a ridge.
	if(tInt>0) ridgePix_A=tInt;
	keeperPix=2*ridgePix_A;
	tInt=round(Dialog.getNumber());//similar, but these are branches from the separated skeleton, so I'll be less stringent
	if(tInt>0) ridgePix_B=tInt;
	//classify ridges
	tReal=Dialog.getNumber();
	if(tReal>=1) maxCurl_line=tReal;
	tReal=Dialog.getNumber();
	if(tReal>=1) minARF_line=tReal;
	tReal=Dialog.getNumber();
	if(tReal>0) dA_ROIend_line=tReal;//I don't think I want to allow a rigorous parallel lines requirement.
	tReal=Dialog.getNumber();
	if(tReal>0) radiusFactor=tReal;
	tReal=Dialog.getNumber();
	if(tReal>0) devLengthRatioThresh=tReal;
	tReal=Dialog.getNumber();
	if(tReal>0) dA_DevAngleDiff=tReal;
	prefReadWrite(false);//causes updated preferences to be written
	if(Dialog.getCheckbox()) printParameters("Log");
}//end doSetupDialog

function printParameters(myDest){
	pString=buildParamString(nl,true);//causes parameter to be on new lines and (true) have labels.
	if(toLowerCase(myDest)=="log") IJ.log(pString);
	else showMessage(pString);
}//end printParameters()

function setUpReport(){
	screenW=screenWidth();
	screenH=screenHeight();
	if(isOpen(liveTableName))print(liveTableNameB,"\\Clear");
	else run("Table...", "name="+liveTableName+" width="+liveTableWide+" height="+liveTableHigh);
	selectWindow(liveTableName);
	setLocation(screenW-(liveTableWide+inset),screenH-liveTableHigh);//shifts table out of the way
	print(liveTableNameB,"\\Headings:"+reportHeader);
}//end setUpReport()

function emptyReport(myID){
	selectImage(myID);
	tTitle=getTitle();
	tStr=tTitle+tab+tab+tab+tab;
	if(!isOpen(liveTableName)) setUpReport();
	print(liveTableNameB,tStr);
}//end emptyReport()

function getROIvalues(myROIindex,myAngles,myRefDist){
	numArrPts=selectionXarr.length;
//ROI location info
	ROI_LocAngle=0;
	ROI_LocDist=0;
	ROI_LocDistRatio=1;
	centerX=getResult("X",myROIindex);
	centerY=getResult("Y",myROIindex);
	if((!isNaN(refX))&&(!isNaN(refY))&&(refX>=0)&&(refY>=0)) {
		dX=centerX-refX;
		dY=centerY-refY;
		ROI_LocAngle=doAngle(dX,dY);
		ROI_LocDist=sqrt(dX*dX+dY*dY);
		if(myRefDist>0) ROI_LocDistRatio=ROI_LocDist/myRefDist;
	}
//ROI  "ends" info
	firstX=selectionXarr[0];
	firstY=selectionYarr[0];
	lastX=selectionXarr[numArrPts-1];
	lastY=selectionYarr[numArrPts-1];
	dX=lastX-firstX;
	dY=lastY-firstY;
	midX=round(firstX+dX/2);
	midY=round(firstY+dY/2);
	endsLength=sqrt(dX*dX+dY*dY);
	endsAngle=doAngle(dX,dY);
	if(abs(endsAngle-ROI_LocAngle)>PI){//adjust endsAngle to be consistent with roi location
		if(ROI_LocAngle<0) endsAngle-=(2*PI);
		else endsAngle+=(2*PI);
	}
//  Profile deviation calculations
//	RMS calculation here based on line connecting end points: http://en.wikipedia.org/wiki/Distance_from_a_point_to_a_line 
//	deviation[ii]=abs(dy*x[ii] - dx*y[ii] + crossProd)/endsLength
	crossProd=lastX*firstY - lastY*firstX;
	sumDev2=0;//sum squares
	maxDev=0;
	maxDevInd=0;
	for(ii=0;ii<numArrPts;ii++){
		tReal=abs(dY*selectionXarr[ii] - dX*selectionYarr[ii] + crossProd)/endsLength;
		if(tReal>maxDev){
			maxDev=tReal;
			maxDevInd=ii;
		}
		sumDev2+=tReal*tReal;
	}
	maxDevRatio=maxDev/endsLength;
	rmsDev=sqrt(sumDev2/numArrPts);
	rmsDevRatio=rmsDev/endsLength;
	devX=selectionXarr[maxDevInd];
	devY=selectionYarr[maxDevInd];
	dX=midX-devX;//set up for vector from bump to bisector of end points
	dY=midY-devY;
	devLength=sqrt(dX*dX+dY*dY);
	devLengthRatio=devLength/endsLength;
	if(devLengthRatio > minDevLengthRatio) devAngle=doAngle(dX,dY);
	else {
		devAngle=endsAngle+PI/2;
		if(devAngle > PI) devAngle -= (2*PI);
		else if(devAngle < PI) devAngle += (2*PI);
	}
	if(abs(devAngle-ROI_LocAngle)>PI){//adjust devAngle to be consistent with roi location
		if(ROI_LocAngle<0) devAngle-=(2*PI);
		else devAngle+=(2*PI);
	}
//Centroid vs ends bisector
	dX=midX-centerX;
	dY=midY-centerY;
	midOffset=sqrt(dX*dX+dY*dY);
	midOffsetRatio=midOffset/endsLength;
//ROI shape calcs
	roiLength=getResult("Perim.",myROIindex)/2;
	roiCurl=roiLength/endsLength;
	roiCurlFeret=roiLength/getResult("Feret",myROIindex);
	AR_feret=getResult("Feret",myROIindex)/getResult("MinFeret",myROIindex);
//UPDATE RESULTS TABLE
	setResult("ROI_LocAngle",myROIindex,ROI_LocAngle);
	setResult("ROI_LocAngleD",myROIindex,ROI_LocAngle*180/PI);
	setResult("ROI_LocDist",myROIindex,ROI_LocDist);
	setResult("ROI_LocDistRatio",myROIindex,ROI_LocDistRatio);
	setResult("EndsLength",myROIindex,endsLength);
	setResult("EndsAngle",myROIindex,endsAngle);
	setResult("EndsAngleD",myROIindex,endsAngle*180/PI);
	setResult("maxDev",myROIindex,maxDev);
	setResult("rmsDev",myROIindex,rmsDev);
	setResult("DevLength",myROIindex,devLength);
	setResult("maxDevRatio",myROIindex,maxDevRatio);
	setResult("rmsDevRatio",myROIindex,rmsDevRatio);
	setResult("DevLengthRatio",myROIindex,devLengthRatio);
	setResult("DevAngle",myROIindex,devAngle);
	setResult("DevAngleD",myROIindex,devAngle*180/PI);
	setResult("MidOffset",myROIindex,midOffset);
	setResult("MidOffsetRatio",myROIindex,midOffsetRatio);
	setResult("ROIlength",myROIindex,roiLength);
	setResult("ROIcurl",myROIindex,roiCurl);
	setResult("ROIcurlF",myROIindex,roiCurlFeret);
	setResult("AR_F",myROIindex, AR_feret);
	setResult("EndsAngleDiff", myROIindex, (endsAngle-ROI_LocAngle)*180/PI);
	setResult("EllipseAngle", myROIindex, getResult("Angle",myROIindex)-ROI_LocAngle*180/PI);
	setResult("DevAngleDiff", myROIindex, (devAngle-ROI_LocAngle)*180/PI);
	return true;
}//end getROIvalues()

function doAngle(dX,dY) {	
	myAngle=NaN;//case of dX=dY=0  should never find this!
	if(dX==0){//check special cases
		if(dY>0) myAngle=-PI/2;
		else if(dY<0) myAngle=PI/2;
		//otherwise dY=0 and angle is undefined
	}
	else if (dY==0){
		if(dX>0) myAngle=0;
		else if(dX<0) myAngle=PI;
	}
	else myAngle=atan2(-dY,dX);//check effect of left-handed coordinate system
	return myAngle;
}//end doAngle()


function getAngleArr(myXarr, myYarr, myPix){//return array of angles at each point in the profile.  
	//myPix is the number of pixels to use around the center pixel:  center pixel +/- (myPix-1)/2
	if(myPix>anglePix) {//rebuild ANGLE_LUT with larger capacity
		if(myPix%2!=1)myPix++;//if not odd, make it odd by making it larger
		anglePix=myPix;//update global
		ANGLE_LUT=makeAngleLUT(anglePix);
	}
	deltaIndexOffset=myPix-1;//n pixels generates delta values up to +/- n-1
	LUTrowLength=2*myPix-1;
	numPts=myXarr.length;
	outArr=newArray(numPts);
	firstIndex=(myPix-1)/2;
	longOffset=firstIndex;
	shortOffset=longOffset-1;
	dXshort=0;
	dYshort=0;
	dXlong=0;
	dYlong=0;
	for(ii=firstIndex;ii<(numPts-firstIndex);ii++){//calculation can only be done on middle of array. Ends handled as special cases
		dXshort=myXarr[ii+shortOffset]-myXarr[ii-shortOffset];
		dYshort=myYarr[ii+shortOffset]-myYarr[ii-shortOffset];
		dXlong=myXarr[ii+longOffset]-myXarr[ii-longOffset];
		dYlong=myYarr[ii+longOffset]-myYarr[ii-longOffset];
		//for three
		xInd=dXshort+deltaIndexOffset;
		yInd=dYshort+deltaIndexOffset;
		LUTthree=yInd*LUTrowLength+xInd;
		//for five
		xInd=dXlong+deltaIndexOffset;
		yInd=dYlong+deltaIndexOffset;
		LUTfive=yInd*LUTrowLength+xInd;
		tThree=ANGLE_LUT[LUTthree];
		tFive=ANGLE_LUT[LUTfive];
		if(abs(tThree-tFive)>PI){//average must "wrap around" +/- PI boundary
			if(tThree<0) tThree+=twoPI;
			else tFive+=twoPI;
			tAngle=(tThree+tFive)/2;
			if(tAngle>PI) tAngle-=twoPI;
			outArr[ii]=tAngle;
		}
		else outArr[ii]=(tThree+tFive)/2;
	}
	fillStart=outArr[firstIndex];
	fillEnd=outArr[numPts-firstIndex-1];
	for(ii=0;ii<firstIndex;ii++) outArr[ii]=fillStart;//fill ends with last known angle
	for(ii=numPts-firstIndex; ii<numPts; ii++) outArr[ii]=fillEnd;
	return outArr;
}//end getAngleArr()

function assessROI(nn) {
	roiType=-1;//-1 for invalid ROI, 0 for line, 1 for chevron
	confFact=0;//essentially a weighting factor where 1 means high confidence about a correct assignment, 0 means no confidence in correct assignment
	lineScore=0;//accumulate score for determining if the blob is a line
	chevronScore=0;//accumulate score for determining if the blob is a line
	midOffsetScore=1;//initialize assuming it will pas
	curlScore=1;//raw end-to-end curl
	AR_FScore=1;//feret based aspect ratio
	//score each factor, add up score:  max = assign roiType=line, confidence = 1
	//first check it for being a line by looking for low curl values.
	currOffsetRatio = getResult("MidOffsetRatio",nn);
	if(currOffsetRatio>midOffsetThresh) midOffsetScore=midOffsetThresh/currOffsetRatio;
	tCurl=getResult("ROIcurl",nn);
	if(tCurl>maxCurl_line) curlScore=maxCurl_line/tCurl;
	tAR_F=getResult("AR_F",nn);
	if(tAR_F<minARF_line)AR_FScore = (tAR_F-1)/minARF_line;
	lineScore=midOffsetScore+curlScore+AR_FScore;
	confFactor=lineScore/lineType_ScoreMax;
	if(lineScore>lineType_ScoreThresh){//check orientation of line
		A_ROIend_score=1;//angle of ends vs location angle of ROI
		tDistRatio=getResult("ROI_LocDistRatio",nn);
		hor_dA_ROIend_line=(1+tDistRatio*radiusFactor)*dA_ROIend_line;
		tROI_LocAngleD=getResult("ROI_LocAngleD",nn);
		quad="UR";//assume 0<=angle<=90
		if(tROI_LocAngleD>90) quad="UL";
		else if(tROI_LocAngleD<=-90) quad="LL";
		else if(tROI_LocAngleD<0) quad="LR";
//		if(quad=="UR"){
			negDiff=-hor_dA_ROIend_line;//start out assuming UR quadrant
			posDiff=dA_ROIend_line;
//		}
//		else if(quad=="UL"){
		if(quad=="UL"){
			negDiff=-dA_ROIend_line;
			posDiff=hor_dA_ROIend_line;
		}
		else if(quad=="LL"){
			negDiff=-hor_dA_ROIend_line;
			posDiff=dA_ROIend_line;
		}
		else if(quad=="LR"){
			negDiff=-dA_ROIend_line;
			posDiff=hor_dA_ROIend_line;
		}
		tAend=getResult("EndsAngleDiff",nn);//this is in degrees
		if((tAend<negDiff)||(tAend>posDiff)) A_ROIend_score = 1-abs(tAend)/90;
		
//		tAend=abs(getResult("EndsAngleDiff",nn));//this is in degrees
//		if(tAend > dA_ROIend_line) A_ROIend_score = 1-tAend/90;
		
		if(A_ROIend_score>=1){
			roiType=1;
			return roiType;
		}
	}
//If not a line, check it for being a chevron pointing correctly
	tDevAngleDiff=abs(getResult("DevAngleDiff",nn));
	tDevLengthRatio=getResult("DevLengthRatio",nn);
	if(tDevLengthRatio > devLengthRatioThresh){//it's gotta look like a chevron!
		if(tDevAngleDiff <= dA_DevAngleDiff) chevronScore=1;
		else chevronScore = 1-tDevAngleDiff/90;
		if(chevronScore > chevron_ScoreThresh){
			roiType = 0;
			return roiType;
		}
	}
//	setResult("EllipseAngle", myROIindex, getResult("Angle",myROIindex)-ROI_LocAngle*180/PI);
//	setResult("DevAngleDiff", myROIindex, (devAngle-ROI_LocAngle)*180/PI);
	return roiType;
}//end assessROI()

function doRidges(myID, myRidgeRad, myRidgeThin, keepRidgeImg) {//returns ridges in ROI Manager, results table is populated with current values
	rtnID=0;
	selectImage(myID);
	run("Duplicate...", "title=tempRidge");
	ridgeID=getImageID();
	run("Find Ridges", "gaussian_radius="+d2s(myRidgeRad,0)+" minimum_thinning="+d2s(myRidgeThin,0));
	setThreshold(2, 255);
	setOption("BlackBackground", false);
	run("Make Binary", "thresholded remaining black");
	run("Options...", "iterations=1 count=1 edm=Overwrite do=Nothing");
	run("Dilate");//heal the "close calls" may need to do this by two iterations
	run("Skeletonize");
	run("Set Measurements...", "area redirect=None decimal=3");
	run("Analyze Particles...", "size="+d2s(ridgePix_A,0)+"-Infinity show=Masks display clear");
	maskID=getImageID();
	run("Invert");
	run("Invert LUT");
	keeperID=isolateKeepers(maskID, keeperPix);//special routine to isolate and keep the "really long" segments.  "Kept" features are put in keeperID and removed from source image
	selectImage(keeperID);
	rename("keepers");
	prunedID=pruneBinary(maskID,ridgePix_B,keepRidgeImg);
	selectImage(keeperID);//now add these to the ROI manager
	run("Set Measurements...", "area redirect=None decimal=3");
	run("Analyze Particles...", "size=0-Infinity show=Nothing display add");//Now should have keepers + pruned in ROImanager
	if(keepRidgeImg) imageCalculator("Min",prunedID,keeperID);//return combined image
	closeByID(maskID);
	closeByID(keeperID);
	closeByID(ridgeID);
	if(!keepRidgeImg) closeByID(prunedID);
	else rtnID=prunedID;
	return rtnID;
}//end doRidges()

function flatten(myID,myOutlierRadius,myOutlierThresh, myMeanRadius){//flattens input image, then smooths
	selectImage(myID);
	run("Remove Outliers...", "radius="+myOutlierRadius+" threshold="+myOutlierThresh+" which=Dark");
	run("Remove Outliers...", "radius="+myOutlierRadius+" threshold="+myOutlierThresh+" which=Bright");
	run("Mean...", "radius="+myMeanRadius);
}//end flatten()

//******************  profile-finding code
function findEnds(myROIindex){//uses roiManager.  Assumes index is valid
	rtnCode=false;
	if(roiManager("count")<=myROIindex) return rtnCode;//check that there is a valid ROI in ROIManager
	roiManager("select",myROIindex);
	run("Interpolate", "interval=1");
	getSelectionCoordinates(xx, yy);
	selectionXarr=xx;//globals, so the calling routine doesn't have to re-extract the values.
	selectionYarr=yy;
	numPts=xx.length-1;//this eliminates the redundant final point in the array
	firstEndIndex=-1;
	secondEndIndex=-1;
	dXfour=0;//initialize these just because
	dYfour=0;
	dXfive=0;
	dYfive=0;
	currX=0; currY=0;
	offset=0;
	for(ii=0;ii<numPts;ii++){//check each point as a candidate for being the beginning of an end!
		currX=xx[ii];
		currY=yy[ii];
		iiPlus2=ii+2;
		if(iiPlus2>=numPts)iiPlus2-=numPts;//wraparound for plus2
		iiPlus4=ii+4;
		if(iiPlus4>=numPts)iiPlus4-=numPts;//wraparound for plus4
		dXfour=xx[iiPlus4]-currX;
		dYfour=yy[iiPlus4]-currY;
		if((dXfour==0)&&(dYfour==0)){//check the easy ones first: lone diagonal pixel, so it's an end
			status=assignEnd(iiPlus2);
			if(status) {
				rtnCode=true;
				break;//break out of ii loop
			}
		}
		else {//inspect the five-step candidates
			iiPlus1=ii+1;
			if(iiPlus1>=numPts)iiPlus1-=numPts;//wraparound for plus1
			iiPlus5=ii+5;
			if(iiPlus5>=numPts)iiPlus5-=numPts;//wraparound for plus5
			dXtwo=xx[iiPlus2]-currX;//need these delta values to satisfy truth table for end patterns
			dYtwo=yy[iiPlus2]-currY;
			dXfive=xx[iiPlus5]-currX;
			dYfive=yy[iiPlus5]-currY;
			if(dXfive==0){//no change in X
				if(dYfive==1){
					if(dXtwo==2){
						if(Roi.contains(xx[ii],yy[ii])) {
							status=assignEnd(iiPlus2);
							if(status) {
								rtnCode=true;
								break;//break out of ii loop
							}
						}
					}
					else if(dXtwo==-2){
						if(Roi.contains(xx[iiPlus1],yy[iiPlus1])) {
							status=assignEnd(iiPlus2);
							if(status) {
								rtnCode=true;
								break;
							}
						}
					}
				}
				else if(dYfive==-1){
					if(dXtwo==2){
						if(Roi.contains(xx[iiPlus4],yy[iiPlus4])) {
							status=assignEnd(iiPlus2);
							if(status) {
								rtnCode=true;
								break;
							}
						}
					}
					else if(dXtwo==-2){
						if(Roi.contains(xx[iiPlus4],yy[iiPlus4])) {
							status=assignEnd(iiPlus2);
							if(status) {
								rtnCode=true;
								break;
							}
						}
					}
				}
			}//end of dXfive=0
			else if(dYfive==0){//no change in Y
				if(dXfive==1){
					if(dYtwo==2){
						if(Roi.contains(xx[ii],yy[ii])) {
							status=assignEnd(iiPlus2);
							if(status) {
								rtnCode=true;
								break;
							}
						}
					}
					else if(dYtwo==-2){
						if(Roi.contains(xx[iiPlus1],yy[iiPlus1])) {
							status=assignEnd(iiPlus2);
							if(status) {
								rtnCode=true;
								break;
							}
						}
					}
				}
				if(dXfive==-1){
					if(dYtwo==2){
						if(Roi.contains(xx[iiPlus4],yy[iiPlus4])) {
							status=assignEnd(iiPlus2);
							if(status) {
								rtnCode=true;
								break;
							}
						}
					}
					else if(dYtwo==-2){
						if(Roi.contains(xx[iiPlus4],yy[iiPlus4])) {
							status=assignEnd(iiPlus2);
							if(status) {
								rtnCode=true;
								break;
							}
						}
					}
				}
			}//end of dXfive=0
		}//end of inspecting the "five step" candidates
	}//end of for() loop
	if(firstEndIndex>secondEndIndex){
		tReal=firstEndIndex;
		firstEndIndex=secondEndIndex;
		secondEndIndex=tReal;
	}
	return rtnCode;
}//end findEnds()

function assignEnd(myIndex){//assumes globals firstEndIndex and secondEndIndex have been initialized
	status=false;
	if(firstEndIndex<0)firstEndIndex=myIndex;
	else {
		secondEndIndex=myIndex;
		status=true;
	}
	return status;
}//end assignEnd()

//*********** ANGLE_LUT **********
function makeAngleLUT(myPix){//returns angles matrix for dX and dY in linear array
	if(myPix%2!=1)myPix++;//if not odd, make it odd by making it larger
	myPts=2*myPix - 1;
	deltaIndexOffset=-(myPts-1)/2;//measured delta can be +/-, so need to offset these so that array indices are in the range 0 to myPts-1
	tArr=newArray(myPts*myPts);
	for(yInd=0;yInd<myPts;yInd++){
		offset=yInd*myPts;
		dY=yInd+deltaIndexOffset;
		for(xInd=0; xInd<myPts; xInd++){
			dX=xInd+deltaIndexOffset;
			index=offset+xInd;
			if((dX==0)&&(dY==0)) tArr[index]=NaN;
			else if(dX==0){
				if(dY>0) tArr[index]=-PI/2;
				else tArr[index]=PI/2;
			}
			else if (dY==0){
				if(dX>0) tArr[index]=0;
				else  tArr[index]=PI;
			}
			else tArr[index]=atan2(-dY,dX);//check effect of left-handed coordinate system
		}
	}
	return tArr;
}//end makeAngleLUT()

function prefReadWrite(myRead){//if myRead==true, then read, otherwise write
	//note: store booleans as "true" or "false", not 1 or 0
	//also note:  order in which they are written matters!!!
	if(myRead){//type arrays are strings
		tStr=call("ij.Prefs.get",STDPREF+keyA, TYPE_A_STR);// pull in Type A, use TYPE_A_STR if item is not found
		result=extractParamString(tStr,sp);//return true if successful, false if it failed
	}
	else {
		aStr=buildParamString(sp,false);//false = no labels
		call("ij.Prefs.set",STDPREF+keyA, aStr);
		call("ij.Prefs.savePreferences");//do this explicitly for safety in case of crash
	}
}//end prefReadWrite

function buildParamString(myDelim,labels){//relies on global variables
	aStr= "FoamPattern_"+DTStamp("yearmonthdayhourminute",uScore);//this is a dummy to make it easier to read the IJ_Prefs
	tStr=d2s(outlierRadius,2);
	if(labels) tStr="outlierRadius"+sp+tStr;
	aStr=aStr+myDelim+tStr;
	tStr=d2s(outlierThresh,0);
	if(labels) tStr="outlierThresh"+sp+tStr;
	aStr=aStr+myDelim+tStr;
	tStr=d2s(meanRadius,2);
	if(labels) tStr="meanRadius"+sp+tStr;
	aStr=aStr+myDelim+tStr;
	tStr="false"; if(useLightRidges) tStr="true";
	if(labels) tStr="useLightRidges"+sp+tStr;
	aStr=aStr+myDelim+tStr;
	tStr="false"; if(useDarkRidges) tStr="true";
	if(labels) tStr="useDarkRidges"+sp+tStr;
	aStr=aStr+myDelim+tStr;
	tStr="false"; if(useVarRidges) tStr="true";
	if(labels) tStr="useVarRidges"+sp+tStr;
	aStr=aStr+myDelim+tStr;
	tStr=d2s(ridgeRad,0);
	if(labels) tStr="ridgeRad"+sp+tStr;
	aStr=aStr+myDelim+tStr;
	tStr=d2s(ridgeThin,0);
	if(labels) tStr="ridgeThin"+sp+tStr;
	aStr=aStr+myDelim+tStr;
	tStr=d2s(varianceRad,2);
	if(labels) tStr="varianceRad"+sp+tStr;
	aStr=aStr+myDelim+tStr;
	tStr=d2s(ridgeRadVar,0);
	if(labels) tStr="ridgeRadVar"+sp+tStr;
	aStr=aStr+myDelim+tStr;
	tStr=d2s(ridgeThinVar,0);
	if(labels) tStr="ridgeThinVar"+sp+tStr;
	aStr=aStr+myDelim+tStr;
	tStr=d2s(ridgePix_A,0);
	if(labels) tStr="ridgePix_A"+sp+tStr;
	aStr=aStr+myDelim+tStr;
	tStr=d2s(ridgePix_B,0);
	if(labels) tStr="ridgePix_B"+sp+tStr;
	aStr=aStr+myDelim+tStr;
	tStr=d2s(maxCurl_line,2);
	if(labels) tStr="maxCurl"+sp+tStr;
	aStr=aStr+myDelim+tStr;
	tStr=d2s(minARF_line,2);
	if(labels) tStr="minARF"+sp+tStr;
	aStr=aStr+myDelim+tStr;
	tStr=d2s(dA_ROIend_line,3);
	if(labels) tStr="dA_ROIend_line"+sp+tStr;
	aStr=aStr+myDelim+tStr;
	tStr=d2s(radiusFactor,3);
	if(labels) tStr="radiusFactor"+sp+tStr;
	aStr=aStr+myDelim+tStr;
	tStr=d2s(devLengthRatioThresh,3);
	if(labels) tStr="devLengthRatioThresh"+sp+tStr;
	aStr=aStr+myDelim+tStr;
	tStr=d2s(dA_DevAngleDiff,3);
	if(labels) tStr="dA_DevAngleDiff"+sp+tStr;
	aStr=aStr+myDelim+tStr;
	tStr=d2s(keeperPix,0);//added keeperPix after this parameter list was already done, so put it at the end.
	if(labels) tStr="keeperPix"+sp+tStr;
	aStr=aStr+myDelim+tStr;
	return aStr;
}//end buildParamString()


function extractParamString(myStr,myDelim){//assumes myStr has no labels
	typeAarr=split(myStr,myDelim);//split string to values,
	test=typeAarr.length;
	ii=0;
	if(ii<test) {tStr=typeAarr[ii]; ii++;}//dummy name
	if(ii<test) {outlierRadius=typeAarr[ii]; ii++;}
	if(ii<test) {outlierThresh=typeAarr[ii]; ii++;}
	if(ii<test) {meanRadius=typeAarr[ii]; ii++;}
	if(ii<test) {if(typeAarr[ii]==trueStr) useLightRidges=true; else useLightRidges=false; ii++;}
	if(ii<test) {if(typeAarr[ii]==trueStr) useDarkRidges=true; else useDarkRidges=false; ii++;}
	if(ii<test) {if(typeAarr[ii]==trueStr) useVarRidges=true; else useVarRidges=false; ii++;}
	if(ii<test) {ridgeRad=typeAarr[ii]; ii++;}
	if(ii<test) {ridgeThin=typeAarr[ii]; ii++;}
	if(ii<test) {varianceRad=typeAarr[ii]; ii++;}
	if(ii<test) {ridgeRadVar=typeAarr[ii]; ii++;}
	if(ii<test) {ridgeThinVar=typeAarr[ii]; ii++;}
	if(ii<test) {ridgePix_A=typeAarr[ii]; ii++;}
	if(ii<test) {ridgePix_B=typeAarr[ii]; ii++;}
	if(ii<test) {maxCurl_line=typeAarr[ii]; ii++;}
	if(ii<test) {minARF_line=typeAarr[ii]; ii++;}
	if(ii<test) {dA_ROIend_line=typeAarr[ii]; ii++;}
	if(ii<test) {radiusFactor=typeAarr[ii]; ii++;}
	if(ii<test) {devLengthRatioThresh=typeAarr[ii]; ii++;}
	if(ii<test) {dA_DevAngleDiff=typeAarr[ii]; ii++;}
	if(ii<test) {keeperPix=typeAarr[ii]; ii++;}	
	return true;
}//end extractParamString

function isolateKeepers(myID, myPix){//warning: resets ROI manager
	/*Special purpose function to keep larger features in a skeletonized image:
		create/remove joints
		check remaining particles for length
			if long enough
				create new image to hold them
				remove them from myID
		return new image containing "keepers"
	*/
	maskID=0;
	selectImage(myID);
	if(!is("binary")) return maskID;//return an impossible ID due to error
	setBackgroundColor(255,255,255);
	setForegroundColor(0,0,0);
	roiManager("Reset")
	roiManager("Show None");
	setOption("BlackBackground", false);
	run("Skeletonize");//should be redundent
	run("Select None");
	run("Duplicate...","title=WorkingImage");
	workID=getImageID();
	run("BinaryConnectivity ", " ");//background is black
	setThreshold(2,3);//select branches and links
	run("Set Measurements...", "area min redirect=None decimal=4");
	run("Analyze Particles...", "size="+d2s(myPix,0)+"-Infinity show=Masks display clear include");
	maskID=getImageID();//this is the image that will be returned as the keepers
	rename("Keepers");
	run("Invert LUT");
	imageCalculator("Max",myID,maskID);
	selectImage(maskID);
	run("Invert");
	closeByID(workID);
	return maskID;
}//end function isolateKeepers()

function pruneBinary(myID, myPrunePix, returnPruned){//returns pruned skeleton image and ID if returnPruned is true.  ROImanager contains pruned skeletons.
	prunedID=0;
	selectImage(myID);
	if(!is("binary")) return prunedID;//return an impossible ID due to error
	setBackgroundColor(255,255,255);
	setForegroundColor(0,0,0);
	roiManager("Reset")
	roiManager("Show None");
	selectImage(myID);
	run("Select None");
	run("Duplicate...","title=WorkingImage");
	initID=getImageID();
	setOption("BlackBackground", false);
	run("Skeletonize");
	linkTitle="LinkMask";
	run("Duplicate...","title="+linkTitle);
	linkID=getImageID();//use this to create the links between nodes
	run("BinaryConnectivity ", " ");//background is black
	nodeTitle="NodeMask";
	run("Duplicate...","title="+nodeTitle);
	nodeID=getImageID();
	setThreshold(4, 9);//select nodes
	run("Make Binary", "thresholded remaining black");
	run("Invert");//nodes are now white on black background
	selectImage(linkID);
	run("Select None");
	setThreshold(2,3);
	run("Set Measurements...", "area min redirect=None decimal=4");
	run("Analyze Particles...", "size=0-Infinity show=Nothing display clear include add");
	numBlobs=roiManager("Count");
	selectImage(nodeID);
	roiManager("Deselect");//make sure no ROIs are selected
	for(ii=numBlobs-1;ii>=0;ii--){
		if(getResult("Min",ii)>2){//if it has no end add it to the nodes (it is part of the backbone)
			roiManager("Select",ii);
			run("Clear");
		}
	}
	imageCalculator("Max", initID,nodeID);
	selectImage(initID);
	run("Select None");
	run("Analyze Particles...", "size="+myPrunePix+"-Infinity show=Masks display include");
	firstPrunedID=getImageID();//contains possible reconstructed lines
	rename("FirstPruned");
	run("Invert");
	run("Invert LUT");
	selectImage(nodeID);
	run("Select None");
	run("Invert");//nodes are now black on white background
	imageCalculator("Min", firstPrunedID,nodeID);
	selectImage(firstPrunedID);
	run("Skeletonize");
	run("BinaryConnectivity ", " ");//background is black
	setThreshold(1, 9);//select all features
	run("Set Measurements...", "area min redirect=None decimal=4");
	run("Analyze Particles...", "size="+myPrunePix+"-Infinity show=Masks display clear include add");//need to filter out junk from isolated nodes
	prunedID=getImageID();//this is ready to move forward
	rename("Pruned");
	run("Invert");
	run("Invert LUT");
	//sort on max bright / remove noded objects
	numNodes=roiManager("Count");
	roiManager("Deselect");//make sure no ROIs are selected
	for(ii=numNodes-1;ii>=0;ii--){
		if(getResult("Max",ii)>3){
			roiManager("Select",ii);
			if(returnPruned)run("Clear");//remove from image only if it is to be returned
			roiManager("Delete");
		}
	}
	selectImage(initID); close();
	selectImage(linkID); close();
	selectImage(nodeID); close();
	selectImage(firstPrunedID); close();
	if(!returnPruned){
		selectImage(prunedID); close();
		prunedID=0;
	}
	return prunedID;
}//end function pruneBinary()

function closeByID(myID){//simplify the check-for-open, then close
	if(isOpen(myID)) {selectImage(myID); close();}
}//end closeByID()

function imageFile(myPath){//tests if file is a valid image extension
	tifExt=".tif"; tiffExt=".tiff";
	jpgExt=".jpg"; jpegExt=".jpeg";
	//what about PNG, BMP, etc.?
	tStr=toLowerCase(myPath);
	if(endsWith(tStr,tifExt)||endsWith(tStr,tiffExt)||endsWith(tStr,jpgExt)||endsWith(tStr,jpegExt))	return true;
	else return false;
}// end imageFile()

function stripExt(inStr,testExt) {//removes extension (or everything after the last ".") from image titles
	lastDotIndex=lastIndexOf(inStr,".");
	if (lastDotIndex>=0){
		rootExt=substring(inStr,lastDotIndex+1,lengthOf(inStr));
		rootFile=substring(inStr,0,lastDotIndex);
	}
	else {
		rootExt="";
		rootFile=inStr;
		return rootFile;
	}
	if(lengthOf(testExt)>0) {//if a test extension was passed, check against it
		if(toLowerCase(rootExt)==toLowerCase(testExt)) return rootFile;
		else return inStr;
	}
	else return rootFile;
}//end stripExt()

function DTStamp(elements,DTSep){
	if(lengthOf(elements)<1)elements="yearmonthdayhourminutesecondmillisecond";//default value if user includes no valid keywords
	else elements = toLowerCase(elements);//convert input to lower case for testing
	elementsList=newArray("year","month","day","hour","minute","second","millisecond");//These are the keys for building result string
	yearInd=0; monthInd=1; dayInd=2; hourInd=3; minInd=4; secInd=5; mSecInd=6;
	yearPos=indexOf(elements,elementsList[yearInd]);
	monthPos=indexOf(elements,elementsList[monthInd]);
	dayPos=indexOf(elements,elementsList[dayInd]);
	hourPos=indexOf(elements,elementsList[hourInd]);
	minPos=indexOf(elements,elementsList[minInd]);
	secPos=indexOf(elements,elementsList[secInd]);
	mSecPos=indexOf(elements,elementsList[mSecInd]);
	orderArray=newArray(yearPos,monthPos,dayPos,hourPos,minPos,secPos,mSecPos);
	Array.sort(orderArray);//now know the final output content and order
	getDateAndTime(year, month, dayOfWeek, dayOfMonth, hour, minute, second, msec);
	datePart="";
	datePartMin=1000;//use these to determine if date part is before/after time part
	timePart="";
	timePartMin=1000;
	gotElement=false;
	for(ii=0;ii<orderArray.length;ii++){
		tInt=orderArray[ii];
		if(tInt >= 0){
			gotElement=true;
			if(tInt==yearPos){datePart=datePart+d2s(year,0); datePartMin=minOf(datePartMin,yearPos);}
			else if (tInt==monthPos){datePart=datePart+IJ.pad(month+1,2); datePartMin=minOf(datePartMin,monthPos);}
			else if (tInt==dayPos){datePart=datePart+IJ.pad(dayOfMonth,2); datePartMin=minOf(datePartMin,dayPos);}
			else if (tInt==hourPos){timePart=timePart+IJ.pad(hour,2); timePartMin=minOf(timePartMin,hourPos);}
			else if (tInt==minPos){timePart=timePart+IJ.pad(minute,2); timePartMin=minOf(timePartMin,minPos);}
			else if (tInt==secPos){timePart=timePart+IJ.pad(second,2); timePartMin=minOf(timePartMin,secPos);}
			else if (tInt==mSecPos){timePart=timePart+IJ.pad(msec,3); timePartMin=minOf(timePartMin,mSecPos);}
		}
	}
	if(!gotElement) return "noDateTimeKeys";
	outString="";
	if(datePartMin<timePartMin)outString=datePart+DTSep+timePart;
	else outString=timePart+DTSep+datePart;
	return outString;
}//end DTStamp()

