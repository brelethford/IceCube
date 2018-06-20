import numpy, pylab, sys, os, healpy

try:
  import cPickle as pickle
except:
  import pickle


import matplotlib


def PlotPValue(ScrambledDir,ObservedTS,Title,Filename,OutDir):
	f=open(ScrambledDir+'/Merged.txt')
	if 'Marginalization' in Title:
		data=2*numpy.log(numpy.array(f.readlines(),dtype=float))
	else:
		data=2*numpy.array(f.readlines(),dtype=float)

	# plot observed TS
	pylab.axvline(ObservedTS,linestyle='solid',linewidth=3,color='red',
		label='Observed test statistic')
	PostTrialsPValue=sum(data>=ObservedTS)/float(len(data))
	pylab.figtext(0.55,0.59,'Post-trial p-value: %5.4f'
		% (round(PostTrialsPValue,3)),fontsize=15)
	#pylab.figtext(0.55,0.59,'Post-trial p-value: %5.4f'
	#	% (round(PostTrialsPValue,2)),fontsize=15)

	# make pretty and save
	pylab.title(r'%s' % (Title),fontsize=20)
	pylab.xlim(xmin=0)
	pylab.ylim(ymax=1.2*max(hist.bincontent))
	pylab.xlabel('Test Statistic',fontsize=18)
	pylab.ylabel('Fraction of trials',fontsize=18)
	pylab.grid(True)
	#pylab.figtext(0.65,0.37,'ICECUBE\nPRELIMINARY',color='red',fontsize=18,
	#	fontweight='bold')
	pylab.legend(prop={'size':'16'})
	pylab.savefig(OutDir+'%s.png' % (Filename))
	pylab.savefig(OutDir+'%s.pdf' % (Filename))
	pylab.clf()
	return PostTrialsPValue

def PlotGalPValue(ScrambledDir,ObservedTS,Title,Filename,OutDir):
	f=open(ScrambledDir+'/Merged.txt')
	if 'Marginalization' in Title:
		data=2*numpy.log(numpy.array(f.readlines(),dtype=float))
	else:
		data=2*numpy.array(f.readlines(),dtype=float)

	# plot observed TS
	pylab.axvline(ObservedTS,linestyle='solid',linewidth=3,color='red',
		label='Observed test statistic')
	PostTrialsPValue=sum(data>=ObservedTS)/float(len(data))
	pylab.figtext(0.60,0.59,'Post-trial p-value: %5.5f'
		% (PostTrialsPValue),fontsize=15)

	# make pretty and save
	pylab.title(r'%s' % (Title),fontsize=20)
	pylab.xlim(xmin=0)
	pylab.xlim(xmax=1.5)
	pylab.ylim(ymax=1.2*max(hist.bincontent))
	pylab.xlabel('Test Statistic',fontsize=18)
	pylab.ylabel('Fraction of trials',fontsize=18)
	pylab.grid(True)
	#pylab.figtext(0.65,0.37,'ICECUBE\nPRELIMINARY',color='red',fontsize=18,
	#	fontweight='bold')
	pylab.legend(prop={'size':'16'})
	pylab.savefig(OutDir+'%s.png' % (Filename))
	#pylab.savefig(OutDir+'%s.eps' % (Filename))
	pylab.clf()
	return PostTrialsPValue

def PlotSourceList(ScrambledDir,Results,Title,Filename,OutDir):
	ObservedTS=2.*max(Results[1])
	PlotPValue(ScrambledDir,ObservedTS,Title,Filename,OutDir)
	# now plot pvalue for each source
	os.system('mkdir -p %s/%s' % (OutDir,Filename))
	PVals=[]
	f=open(OutDir+Filename+'.txt','w')
	f.write('SourceName & TestStatistic & BestFitNs & Pre-trials P-value\n')
	i=0
	while i < len(Results[0]):
		source=Results[0][i]
		ObservedTS=2.*Results[1][i]
		PVals+=[PlotPValue(ScrambledDir+'/%s/' % (source),ObservedTS,source,
			source,OutDir+'/'+Filename+'/')]
		f.write('%s & %5.4f & %5.2f & %5.3f\n' % (Results[0][i],
			2.*Results[1][i],Results[2][i],PVals[i]))
		i+=1
	ind = numpy.argmax(Results[1])
	f.write('\n')
	f.write('Source with highest TS: %s %5.4f %5.4f\n' % (Results[0][ind],
		2.*Results[1][ind],Results[2][ind]))
	f.close()

def PlotGalacticPlaneScan(ScrambledDir,ResultDict,Filename,OutDir):
	# first plot TS plots for each width
	os.system('mkdir -p %s/GalacticPlaneScan' % (OutDir))
	ScrambledTS=[]
	ObservedTS=[]
	Widths=[]
	PVals=[]
	for item in ResultDict:
		width=item.split('_')[-1]
		res=ResultDict[item]
		TS=2*res[0]
		title='Galactic Plane, Width = %s$^{\circ}$' % (width)
		pval=PlotPValue(ScrambledDir+'/'+width,TS,title,
			Filename+'_%s' % (width),OutDir+'/GalacticPlaneScan/')
		Widths+=[float(width)]
		PVals+=[pval]
		ObservedTS+=[TS]
		ScrambledTS+=[numpy.genfromtxt(ScrambledDir+'/'+width+'/Merged.txt')]

	# calculate and plot post-trials p-value
	ScrambledTS=numpy.transpose(ScrambledTS)
	ScrambledTSMax=numpy.max(ScrambledTS,axis=1)
	assert len(ScrambledTSMax)==len(ScrambledTS[:,0])
	f=open(ScrambledDir+'/Merged.txt','w')
	for item in ScrambledTSMax:
		f.write('%3.4f\n' % (item))
	f.close()
	title='Galactic Plane Scan, Post-trials Test Statistic'
	PostTrialsPValue=PlotPValue(ScrambledDir,max(ObservedTS),title,Filename,OutDir)

	assert PostTrialsPValue==1.*sum(2*ScrambledTSMax>max(ObservedTS))/\
		len(ScrambledTSMax)

	# now plot local p-value vs. width
	Widths=numpy.array(Widths)
	PVals=numpy.array(PVals)
	pylab.scatter(Widths,PVals,facecolor='black')
	pylab.xlim(0,30)
	pylab.ylim(0,max(PVals))
	title=pylab.title(r'Galactic Plane with $|b| < \theta_{max}$',fontsize=20)
	title.set_y(1.01)
	pylab.xlabel(r'$\theta_{max}[^{\circ}]$',fontsize=18)
	#pylab.xlabel(r'Width of Galactic Plane Hypothesis ($^{\circ}$)')
	pylab.ylabel('Pre-trials p-value',fontsize=18)
	pylab.savefig(OutDir+Filename+'_PVal_Vs_Width_NoFigText.pdf')
	pylab.savefig(OutDir+Filename+'_PVal_Vs_Width_NoFigText.png')
	pylab.figtext(0.525,0.75,'Post-trials p-value: %5.4f'
		% (PostTrialsPValue),fontsize=17)
	pylab.figtext(0.2,0.75,'ICECUBE\nPRELIMINARY',color='red',fontsize=16,
		fontweight='bold')
	pylab.savefig(OutDir+Filename+'_PVal_Vs_Width.pdf')
	pylab.savefig(OutDir+Filename+'_PVal_Vs_Width.png')
	pylab.clf()

def AddEventMarkersAndNumbers(cascadesonly=False, energythres=0.0,mycoord=['G'],Galactic=False):
	
	# plot and label event locations using unsmoothed positions from claudio
	# TODO: edit this for new events
	#	-make sure labels all appear in a good place
	
	#load data
	datadir = '/data/user/wgiang/yr4HESE/studies/'
	datafile = 'eventinfo_fixed.txt'

	allIDs, allenergies, allDECs, allRAs, alleventtypes = numpy.loadtxt(
                datadir+datafile,
                dtype={'names':('ID','Energy','dec','ra','eventtype'),
                       'formats':('f2','f7','f5','f5','S8')},
                usecols=(0,1,3,4,5),
                unpack=True)

	assert(len(allenergies) == len(allRAs) == len(allDECs) == 52) #4yrHESE

	cut = (allenergies >= energythres) #TeV
	if cascadesonly:
		cut = cut & (alleventtypes == 'Shower')

	energies = allenergies[cut]
	RAs = allRAs[cut]
	DECs = allDECs[cut]
	eventtypes = alleventtypes[cut]
	tracks = eventtypes == 'Track'
	cascades = eventtypes == 'Shower'
	IDs = allIDs[cut]
	
	thetas = numpy.radians(DECs-90.0)
	phis = numpy.radians(RAs+180)%(2*numpy.pi)

	# after removing events 28 and 32, this is how the skymaps are ordered.
	name = numpy.array(['3','12','14','10','25','7','16','5','20','8','1','23','21','26',\
		'17','19','22','13','9','6','24','11','27','2','4','15','18','29',\
		'30','31','33','34','35','36','37','38','39','40','41','42', '43',\
		'44', '45','46','47','48','49','50','51','52','53','54'])
	
	IDs_with_cut = numpy.array([str(int(ID)) for ID in IDs])

	passescut = [numpy.where(name == ID)[0][0] for ID in IDs_with_cut]
	name = name[passescut] # Now name is ordered to match DECs/RAs
	
	offsetE = numpy.zeros(sum(cascades))
	offsetE[3] = 0.02
	if sum(tracks>0):
		healpy.projscatter(thetas[tracks], phis[tracks],
			s=50,marker='x', color='k',coord=mycoord)
	healpy.projscatter(thetas[cascades],phis[cascades]+offsetE,s=90,
		marker='+',color='k',coord=mycoord)

	if Galactic:
		# mark north/south
		healpy.projtext(64,-17,s='Northern Hemisphere',coord=['G'],lonlat=True,
			fontsize=7.5,rotation=71)
		healpy.projtext(57,-19,s='Southern Hemisphere',coord=['G'],lonlat=True,
			fontsize=7.5,rotation=71)
		# use healpy to rotate event coords from equatorial to galactic, then
		# add offsets to make the plots pretty
		transform=healpy.Rotator(coord=['C','G'])
		phig,thetag = transform(RAs,DECs,lonlat=True)
		for i in xrange(len(name)):
			if name[i] in ['13','33','25','12']:
				healpy.projtext(phig[i]+12, thetag[i]-2, fontsize=11,
					s=name[i],coord=['G'],lonlat=True)
			elif name[i] in ['26','9']:
				healpy.projtext(phig[i]+18, thetag[i]-2, fontsize=11,
					s=name[i],coord=['G'],lonlat=True)
			elif name[i] in ['37']:
				healpy.projtext(phig[i]+25, thetag[i]-10, fontsize=11,
					s=name[i],coord=['G'],lonlat=True)
			elif name[i] in ['3']:
				healpy.projtext(phig[i]+9, thetag[i]-2, fontsize=11,
					s=name[i],coord=['G'],lonlat=True)
			elif name[i] in ['14']:
				healpy.projtext(phig[i]-4, thetag[i]+2, fontsize=11,
					s=name[i],coord=['G'],lonlat=True)
			elif name[i] in ['18']:
				healpy.projtext(phig[i]+15, thetag[i]+5, fontsize=11,
					s=name[i],coord=['G'],lonlat=True)
			elif name[i] in ['21']:
				healpy.projtext(phig[i]+9, thetag[i]+5, fontsize=11,
					s=name[i],coord=['G'],lonlat=True)
			elif name[i] in ['10']:
				healpy.projtext(phig[i]+12, thetag[i]+7, fontsize=11,
					s=name[i],coord=['G'],lonlat=True)
			elif name[i] in ['8']:
				healpy.projtext(phig[i]-3, thetag[i]+3, fontsize=11,
					s=name[i],coord=['G'],lonlat=True)
			elif name[i] in ['38','44']:
				healpy.projtext(phig[i]-3, thetag[i]-4, fontsize=11,
					s=name[i],coord=['G'],lonlat=True)
			elif name[i] in ['48']:
				healpy.projtext(phig[i]+12, thetag[i]-2,fontsize=11,
					s=name[i],coord=['G'],lonlat=True)
			elif name[i] in ['49']:
				healpy.projtext(phig[i]-1, thetag[i]-4, fontsize=11,
					s=name[i],coord=['G'],lonlat=True)
			else:
				healpy.projtext(phig[i]+3, thetag[i]+3, fontsize=11,
					s=name[i],coord=['G'],lonlat=True)
	else:
		thetag = thetas
		phig   = phis

		for i in range(len(name)):
			if name[i]=='14' or name[i]=='21'\
				or name[i]=='10' or name[i]=='22':
				healpy.projtext(thetag[i], phig[i]+0.2, fontsize=11,
					s=name[i],coord=mycoord)
			elif name[i]=='2':
				healpy.projtext(thetag[i], phig[i]+0.1, fontsize=11,
					s=name[i],coord=mycoord)
			elif name[i]=='12':
				healpy.projtext(thetag[i], phig[i]+0.25, fontsize=11,
					s=name[i], coord=mycoord)
			elif name[i]=='13' or name[i]=='17' or name[i]=='26':
				healpy.projtext(thetag[i]+0.03, phig[i]-0.05,
					fontsize=11,s=name[i], coord=mycoord)
			#elif name[i]=='23':
				#healpy.projtext(thetag[i]-0.1, phig[i]+0.16,
					#fontsize=11,s=name[i], coord=mycoord)
			elif name[i]=='3':
				healpy.projtext(thetag[i]-0.1, phig[i]+0.1,
					fontsize=11,s=name[i], coord=mycoord)
			elif name[i]=='5':
				healpy.projtext(thetag[i]+0.05, phig[i]+0.05,
					fontsize=11,s=name[i], coord=mycoord)
			elif name[i]=='15':
				healpy.projtext(thetag[i]-0.05, phig[i]-0.02,
					fontsize=11,s=name[i], coord=mycoord)
			elif name[i]=='24':
				healpy.projtext(thetag[i]+0.04, phig[i]-0.04,
					fontsize=11,s=name[i], coord=mycoord)
			elif name[i]=='25':
				healpy.projtext(thetag[i]+0.03, phig[i]+0.15,
					fontsize=11,s=name[i], coord=mycoord)
			elif name[i]=='31':
				healpy.projtext(thetag[i]-0.05, phig[i]-0.1,
					fontsize=11,s=name[i], coord=mycoord)
			elif name[i]=='29' or name[i]=='34':
				healpy.projtext(thetag[i]-0.05, phig[i]-0.13,
					fontsize=11,s=name[i], coord=mycoord)
			elif name[i]=='33':
				healpy.projtext(thetag[i]+0.05, phig[i]+0.1,
					fontsize=11,s=name[i], coord=mycoord)
			elif name[i]=='35':
				healpy.projtext(thetag[i]+0.1, phig[i]+0.1,
					fontsize=11,s=name[i], coord=mycoord)
			elif name[i]=='36':
				healpy.projtext(thetag[i]+0.05, phig[i]+0.1,
					fontsize=11,s=name[i], coord=mycoord)
			elif name[i]=='37':
				healpy.projtext(thetag[i]+0.05, phig[i]+0.15,
					fontsize=11,s=name[i], coord=mycoord)
			elif name[i]=='16':
				healpy.projtext(thetag[i]-0.12, phig[i]+0.02,
					fontsize=11,s=name[i], coord=mycoord)
			elif name[i]=='8':
				healpy.projtext(thetag[i]+0.03, phig[i]-0.04,
					fontsize=11,s=name[i], coord=mycoord)
			elif name[i]=='43':
				healpy.projtext(thetag[i]+0.04, phig[i]-0.02,
					fontsize=11,s=name[i], coord=mycoord)
			elif name[i]=='48' or name[i]=='49':
				healpy.projtext(thetag[i]-0.12, phig[i]+0.04,
					fontsize=11,s=name[i], coord=mycoord)
			else:
				healpy.projtext(thetag[i]+0.05, phig[i]+0.1,
					fontsize=11,s=name[i], coord=mycoord)

def PlotSkyMap(LLHMap,Filename,OutDir,PlotSourceList=False,Galactic=False,
	NoPreliminary=False,iscascadesonly=False,energythreshold=0.0):
	# color scheme to look pretty
	scale = 0.9999
	offset = 1-scale
	mydict = {'red': ((0.0, 1.0, 1.0),
		 (0.0*scale+offset, 0.9961, 0.9961),
		 (0.5*scale+offset, 0.9235, 0.9235),
		 (1.0*scale+offset, 0.6218, 0.6218),),
	'green': ((0.0, 1.0, 1.0),
		(0.0*scale+offset, 0.9725, 0.9725),
		(0.5*scale+offset, 0.5000, 0.5000),
		(1.0*scale+offset, 0.1753, 0.1753),),
	'blue': ((0.0, 1.0, 1.0),
		(0.0*scale+offset, 0.9961, 0.9961),
		(0.5*scale+offset, 0.9341, 0.9341),
		(1.0*scale+offset, 1.0000, 1.0000),)}
	mycm = matplotlib.colors.LinearSegmentedColormap('my_ncar',mydict, 256)

	# plot skymap
	if Galactic:
		healpy.mollview(2*LLHMap[::-1],title="", coord=['C','G'], cmap=mycm,
			min=0.0,max=round(max(2*LLHMap),1),unit='TS=2log(L/L0)')
		healpy.graticule(coord=['G'],color='DimGrey')
		mycoord=['C','G']
		healpy.projtext(numpy.pi/2-0.025, -numpy.pi+0.4, s='-180$^\circ$',
			color='DimGrey',coord=['G'])
		healpy.projtext(numpy.pi/2-0.025, numpy.pi-0.05, color='DimGrey',
			coord=['G'], s='180$^\circ$')
		# plot equatorial plane
		horizon = 90.
		ras = numpy.arange(0.,361.,1.)*numpy.pi/180.
		decls_1 = numpy.ones(len(ras))*(180.-horizon) *numpy.pi/180.
		healpy.projplot(decls_1,ras, color='DimGrey', linewidth=1., alpha=0.5,
			coord=['C'],)
		decls_2 = numpy.ones(len(ras))*(180.-horizon+2.5) *numpy.pi/180.
		decls_3 = numpy.ones(len(ras))*(180.-horizon-2.5) *numpy.pi/180.

		healpy.projplot(decls_2,ras,color='red',linewidth=1,coord=['G'])
	
		healpy.projplot(decls_3,ras,color='red',linewidth=1,coord=['G'])
	else:
		healpy.mollview(2*LLHMap[::-1],title="", coord=['C'], cmap=mycm,
			min=0.0,max=round(max(2*LLHMap),1),unit='TS=2log(L/L0)', rot=180)
		healpy.graticule(coord=['C'],color='DimGrey')
		healpy.projtext(numpy.pi/2, 0.15, s='0$^\circ$', color='DimGrey')
		healpy.projtext(numpy.pi/2, 2*numpy.pi-0.05, color='DimGrey',
			s='360$^\circ$')
		mycoord=['C']

		# plot galactic plane and center
		horizon = 90.
		ras = numpy.arange(0.,361.,1.)*numpy.pi/180.
		decls_1 = numpy.ones(len(ras))*(180.-horizon) *numpy.pi/180.
		healpy.projplot(decls_1,ras, color='DimGrey', linewidth=1., alpha=0.5,
			coord=['G'],)
		healpy.projscatter(0.0, 0.0, color='DimGrey', marker='s',coord=['G'],
			lonlat=True)
	
	if PlotSourceList:
		sl=numpy.genfromtxt('SL_All.txt',dtype=None)
		healpy.projscatter(sl['f1'],sl['f2'],lonlat=True,s=20,marker='o',
			color='black',coord=mycoord)

	AddEventMarkersAndNumbers(iscascadesonly, energythreshold, mycoord=mycoord,
		Galactic=Galactic)

	if not NoPreliminary:
		pylab.figtext(0.05,0.925,'ICECUBE PRELIMINARY',color='red',fontsize=14,fontweight='bold')

	# save
	pylab.savefig(OutDir+Filename+'.png')
	pylab.savefig(OutDir+Filename+'.pdf')
	pylab.clf()


