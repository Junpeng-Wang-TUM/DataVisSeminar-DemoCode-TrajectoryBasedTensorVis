%% This code is created for visualizing the 2D stress tensor field via the evenly-spaced Principal Stress Lines
%% Author: Junpeng Wang (junpeng.wang@tum.de)
%% Date: 2021-09-14
%% Refer to Fig.2b in paper "Stress Topology Analysis for Porous Infill Optimization" (2021, arXiv:2108.09675) by 
%%	Wang, J., Wu, J. and Westermann, R. for introductory understanding.
%% The method used here is from the submission paper "The 3D Trajectory-based Stress Visualizer" (2021, under review) by
%%	Junpeng Wang, Christoph Neuhauser, Jun Wu, Xifeng Gao and Rüdiger Westermann,

function TSV2D(stressfileName,PSLsDensityCtrl)
	global snappingOpt_;
	
	%%1. Read Data
	LoadStressField(stressfileName);
	
	%%2. Setup
	Setup();

	%%3. Topology Analysis or PSLs Generation
	figure(1); 
	snappingOpt_ = 0;
	CreatePSLsViaTSV(PSLsDensityCtrl);
	ShowPSLs();
end

function Setup()
	global eleSize_;
	global tracingStepWidth_;
	global relaxedFactor_;
	global degeneratePoints_;
	global majorPSLpool_;
	global minorPSLpool_;
	
	relaxedFactor_ = 1.0;
	majorPSLpool_ = PrincipalStressLineStruct();
	minorPSLpool_ = PrincipalStressLineStruct();	
	tracingStepWidth_ = eleSize_;
end

function LoadStressField(fileName)
	global boundingBox_;
	global numNodes_;
	global nodeCoords_;
	global numEles_;
	global eNodMat_;
	global eleState_;
	global nodState_;
	global elementsOnBoundary_;
	global nodesOnBoundary_;
	global cartesianStressField_;
	global loadingCond_; 
	global fixingCond_;
	global eleCentroidList_;

	global nelx_;
	global nely_;
	global carNodMapForward_;
	global voxelizedVolume_;	
	
	%%read mesh and cartesian stress field
	fid = fopen(fileName, 'r');
	fgetl(fid); fgetl(fid); fgetl(fid); 
	tmp = fscanf(fid, '%s', 1);	
	meshType = fscanf(fid, '%s', 1);
	
	tmp = fscanf(fid, '%s', 1);
	tmp = fscanf(fid, '%d %d', [1 2]);
	nelx_ = tmp(1); nely_ = tmp(2); 
	tmp = fscanf(fid, '%s', 1);
	lBound = fscanf(fid, '%f %f', [1 2]);
	tmp = fscanf(fid, '%s', 1);
	uBound = fscanf(fid, '%f %f', [1 2]);
	boundingBox_ = [lBound; uBound];
	tmp = fscanf(fid, '%s', 1); 
	numValidEles = fscanf(fid, '%d', 1);
	tmp = fscanf(fid, '%s', 1);
	validElements = fscanf(fid, '%d', [1, numValidEles])';
	validElements = validElements + 1;		

	%%read cartesian stress field
	tmp = fscanf(fid, '%s %s %s %s', 4);
	numStressFields = fscanf(fid, '%d', 1);
	tmp = fscanf(fid, '%s %s', 2); numLoadedNodes = fscanf(fid, '%d', 1);
	tmp = fscanf(fid, '%d %f %f', [3, numLoadedNodes]); 
	tmp(1,:) = tmp(1,:)+1; 
	loadingCond_ = tmp';
	tmp = fscanf(fid, '%s %s', 2); numFixedNodes = fscanf(fid, '%d', 1);
	tmp = fscanf(fid, '%d', [1, numFixedNodes]); 
	fixingCond_ = tmp'+1;
	tmp = fscanf(fid, '%s %s', 2); numValidNods = fscanf(fid, '%d', 1);
	tmp = fscanf(fid, '%f %f %f', [3, numValidNods]);
	cartesianStressField_ = tmp';
	fclose(fid);

	%%recover cartesian mesh
	voxelizedVolume_ = zeros(nelx_*nely_,1);
	voxelizedVolume_(validElements) = 1;
	voxelizedVolume_ = reshape(voxelizedVolume_, nely_, nelx_);
	RecoverCartesianMesh();
	numNod2ElesVec = zeros(numNodes_,1);
	for ii=1:numEles_
		iNodes = eNodMat_(ii,:);
		numNod2ElesVec(iNodes,1) = numNod2ElesVec(iNodes,1)+1;
	end
	nodesOnBoundary_ = find(numNod2ElesVec<4);	
	nodState_ = zeros(numNodes_,1); nodState_(nodesOnBoundary_) = 1;
	eleState_ = 4*ones(1, numEles_);	
	%%%%
	allNodes = zeros(numNodes_,1);
	allNodes(nodesOnBoundary_) = 1;	
	tmp = zeros(numEles_,1);
	for ii=1:4
		tmp = tmp + allNodes(eNodMat_(:,ii));
	end
	elementsOnBoundary_ = find(tmp>0);				
	%%%%
	%% element centroids
	eleNodCoordListX = nodeCoords_(:,1); eleNodCoordListX = eleNodCoordListX(eNodMat_);
	eleNodCoordListY = nodeCoords_(:,2); eleNodCoordListY = eleNodCoordListY(eNodMat_);
	eleCentroidList_ = [sum(eleNodCoordListX,2) sum(eleNodCoordListY,2)]/4;
	
	loadingCond_(:,1) = carNodMapForward_(loadingCond_(:,1));
	fixingCond_ = carNodMapForward_(fixingCond_);
end

function RecoverCartesianMesh()	
	global nelx_; global nely_;
	global voxelizedVolume_;
	global boundingBox_;
	global numEles_; global numNodes_; global eleSize_;
	global nodeCoords_; global eNodMat_; 
	global carEleMapBack_; global carEleMapForward_;
	global carNodMapBack_; global carNodMapForward_;
	global meshState_; 
	%    __ x
	%   / 
	%  -y         
	%		4--------3
	%	    |		 |		
	%		|		 |
	%		1--------2
	%	rectangular element		
	eleSize_ = min((boundingBox_(2,:)-boundingBox_(1,:))./[nelx_ nely_]);
	carEleMapBack_ = find(1==voxelizedVolume_);	
	numEles_ = length(carEleMapBack_);
	meshState_ = zeros(nelx_*nely_,1);	
	meshState_(carEleMapBack_) = 1;	
	carEleMapForward_ = zeros(nelx_*nely_,1);	
	carEleMapForward_(carEleMapBack_) = (1:numEles_)';
	nodenrs = reshape(1:(nelx_+1)*(nely_+1), 1+nely_, 1+nelx_);
	eNodVec = reshape(nodenrs(1:end-1,1:end-1)+1, nelx_*nely_, 1);
	eNodMat_ = repmat(eNodVec(carEleMapBack_),1,4);
	tmp = [0 nely_+[1 0] -1];
	for ii=1:4
		eNodMat_(:,ii) = eNodMat_(:,ii) + repmat(tmp(ii), numEles_,1);
	end
	carNodMapBack_ = unique(eNodMat_);
	numNodes_ = length(carNodMapBack_);
	carNodMapForward_ = zeros((nelx_+1)*(nely_+1),1);
	carNodMapForward_(carNodMapBack_) = (1:numNodes_)';	
	for ii=1:4
		eNodMat_(:,ii) = carNodMapForward_(eNodMat_(:,ii));
	end
	nodeCoords_ = zeros((nelx_+1)*(nely_+1),2);

	xSeed = boundingBox_(1,1):(boundingBox_(2,1)-boundingBox_(1,1))/nelx_:boundingBox_(2,1);
	ySeed = boundingBox_(2,2):(boundingBox_(1,2)-boundingBox_(2,2))/nely_:boundingBox_(1,2);		
	nodeCoords_(:,1) = reshape(repmat(xSeed, nely_+1, 1), (nelx_+1)*(nely_+1), 1);
	nodeCoords_(:,2) = repmat(ySeed, 1, nelx_+1)';	
	
	nodeCoords_ = nodeCoords_(carNodMapBack_,:);
end

function N = ShapeFunction(s, t)
	%				   	   __s (parametric coordinate system)
	%				  	 /-t
	%				*4			*3
	%			*1			*2
	%
	%				nodes
	s = s(:);
	t = t(:);
	N = zeros(size(s,1), 4);
	N(:,1) = 0.25*(1-s).*(1-t);
	N(:,2) = 0.25*(1+s).*(1-t);
	N(:,3) = 0.25*(1+s).*(1+t);
	N(:,4) = 0.25*(1-s).*(1+t);	
end

function val = PrincipalStressLineStruct()
	val = struct(...
		'ith',						0, 	...
		'length',					0,	...
		'midPointPosition',			0,	...		
		'phyCoordList',				[], ...
		'eleIndexList',				[], ...
		'principalStressList',		[] ...
	);	
end

function principalStress = ComputePrincipalStress(cartesianStress)
	principalStress = zeros(size(cartesianStress,1), 1+2+1+2);
	iPS = zeros(1, 6);
	for ii=1:size(cartesianStress,1)
		iCartesianStress = cartesianStress(ii,:);
		A = iCartesianStress([1 3; 3 2]);
		[eigenVec, eigenVal] = eig(A);
		iPS([1 4]) = diag(eigenVal);
		iPS([2 3 5 6]) = reshape(eigenVec,1,4);
		principalStress(ii,:) = iPS;
	end		
end

function [majorPSL, minorPSL] = GeneratePrincipalStressLines(initialSeed, iniDir, limiSteps, typePSL)
	majorPSL = PrincipalStressLineStruct();
	minorPSL = PrincipalStressLineStruct();
	
	%%1. Spot the Starting Point	
	[eleIndex, phyCoord, principalStress, opt] = PreparingForTracing(initialSeed);
	if ~opt, return; end
	%%2. Compute PSL(s)
	if ~isempty(iniDir)
		psDir = [5 6];
		majorPSL = ComputePSL(eleIndex, phyCoord, principalStress, psDir, iniDir, limiSteps);
		psDir = [2 3];
		minorPSL = ComputePSL(eleIndex, phyCoord, principalStress, psDir, iniDir, limiSteps);
	else
		switch typePSL
			case 'MAJOR'
				psDir = [5 6];
				majorPSL = ComputePSL(eleIndex, phyCoord, principalStress, psDir, principalStress(psDir), limiSteps);			
			case 'MINOR'
				psDir = [2 3];			
				minorPSL = ComputePSL(eleIndex, phyCoord, principalStress, psDir, principalStress(psDir), limiSteps);			
			case 'BOTH'
				psDir = [5 6];
				majorPSL = ComputePSL(eleIndex, phyCoord, principalStress, psDir, principalStress(psDir), limiSteps);	
				psDir = [2 3];
				minorPSL = ComputePSL(eleIndex, phyCoord, principalStress, psDir, principalStress(psDir), limiSteps);					
		end
	end
end

function iPSL = ComputePSL(eleIndex, phyCoord, principalStress, psDir, iniDir, limiSteps)
	global tracingStepWidth_;
	iPSL = PrincipalStressLineStruct();
	PSLphyCoordList = phyCoord;
	PSLeleIndexList = eleIndex;
	PSLprincipalStressList = principalStress;
	%% tracing along first direction (v1)
	nextPoint = phyCoord + tracingStepWidth_*iniDir;
	[phyCoordList, eleIndexList, principalStressList] = TracingPSL(nextPoint, iniDir, eleIndex, psDir, limiSteps);
	PSLphyCoordList = [PSLphyCoordList; phyCoordList];
	PSLeleIndexList = [PSLeleIndexList; eleIndexList];
	PSLprincipalStressList = [PSLprincipalStressList; principalStressList];
	%% tracing along second direction (-v1)			
	nextPoint = phyCoord - tracingStepWidth_*iniDir;
	[phyCoordList, eleIndexList, principalStressList] = TracingPSL(nextPoint, -iniDir, eleIndex, psDir, limiSteps);
	if size(phyCoordList,1) > 1
		phyCoordList = flip(phyCoordList);
		eleIndexList = flip(eleIndexList);
		principalStressList = flip(principalStressList);				
	end
	PSLphyCoordList = [phyCoordList; PSLphyCoordList];
	PSLeleIndexList = [eleIndexList; PSLeleIndexList];
	PSLprincipalStressList = [principalStressList; PSLprincipalStressList];
	iPSL.midPointPosition = size(phyCoordList,1)+1;	
	%%2.3 finish Tracing the current major PSL			
	iPSL.length = size(PSLphyCoordList,1);
	iPSL.eleIndexList = PSLeleIndexList;
	iPSL.phyCoordList = PSLphyCoordList;
	iPSL.principalStressList = PSLprincipalStressList;				
end

function [eleIndex, phyCoord, principalStress, opt] = PreparingForTracing(initialSeed)
	global nodeCoords_; global eNodMat_;
	global cartesianStressField_;
	eleIndex = 0;
	phyCoord = 0; 
	principalStress = 0;
	if 3==size(initialSeed,2)
		opt = 1;
		formatedSeed = initialSeed;
		eleIndex = formatedSeed(1,1);
	elseif 2==size(initialSeed,2)
		[eleIndex, paraCoordinates, opt] = FindAdjacentElement(initialSeed);		
		if opt
			formatedSeed = [eleIndex paraCoordinates];
		else
			return;
		end
	else
		error('Wrong Input!');
	end
	NIdx = eNodMat_(eleIndex,:)';
	eleNodeCoords = nodeCoords_(NIdx,:);
	eleCartesianStress = cartesianStressField_(NIdx,:);
	paraCoord = formatedSeed(1, 2:3);
	shapeFuncs = ShapeFunction(paraCoord(1), paraCoord(2));	
	phyCoord = shapeFuncs*eleNodeCoords;						
	interpolatedCartesianStress = shapeFuncs*eleCartesianStress;
	principalStress = ComputePrincipalStress(interpolatedCartesianStress);	
end

function [phyCoordList, eleIndexList, principalStressList] = TracingPSL(nextPoint, iniDir, elementIndex, typePSL, limiSteps)
	%% Tracing the PSL by 2-nd order Runge-Kutta Scheme 
	global eNodMat_;
	global cartesianStressField_;
	global tracingStepWidth_; 
	
	phyCoordList = zeros(limiSteps,2);
	eleIndexList = zeros(limiSteps,1);
	principalStressList = zeros(limiSteps,6);
	index = 0;	
	
	intgerScheme = 'RK2'; %% 'RK2', 'EULER'
	switch intgerScheme		
		case 'EULER'
			[elementIndex, paraCoordinates, bool1] = FindAdjacentElement(nextPoint);		
			while 1==bool1
				index = index + 1; if index > limiSteps, index = index-1; break; end
				cartesianStress = cartesianStressField_(eNodMat_(elementIndex,:)', :);
				shapeFuncs = ShapeFunction(paraCoordinates(1), paraCoordinates(2));
				cartesianStressOnGivenPoint = shapeFuncs*cartesianStress;
				principalStress = ComputePrincipalStress(cartesianStressOnGivenPoint);					
				nextDir = DirectionSelecting(iniDir, principalStress(typePSL), -principalStress(typePSL));
					
				if 0 == AngleTerminationCondition(iniDir, nextDir), index = index-1; break; end	
				iniDir = nextDir;
				phyCoordList(index,:) = nextPoint;
				eleIndexList(index,:) = elementIndex;
				principalStressList(index,:) = principalStress;					
				nextPoint = nextPoint + tracingStepWidth_*iniDir;
				[elementIndex, paraCoordinates, bool1] = FindAdjacentElement(nextPoint);
			end				
		case 'RK2'
			%%initialize initial k1 and k2
			k1 = iniDir;
			iniPot = nextPoint - k1*tracingStepWidth_;
			midPot = nextPoint - k1*tracingStepWidth_/2;
			
				
			[elementIndex, paraCoordinates, bool1] = FindAdjacentElement(midPot);
			if bool1
				cartesianStress = cartesianStressField_(eNodMat_(elementIndex,:)', :);
				shapeFuncs = ShapeFunction(paraCoordinates(1), paraCoordinates(2));
				cartesianStressOnGivenPoint = shapeFuncs*cartesianStress;
				principalStress = ComputePrincipalStress(cartesianStressOnGivenPoint);
				k2 = DirectionSelecting(k1, principalStress(typePSL), -principalStress(typePSL));
				nextPoint = iniPot + tracingStepWidth_*k2;
				[elementIndex, paraCoordinates, bool1] = FindAdjacentElement(nextPoint);
				while 1==bool1
					index = index + 1; if index > limiSteps, index = index-1; break; end
					%%k1
					cartesianStress = cartesianStressField_(eNodMat_(elementIndex,:)', :);
					shapeFuncs = ShapeFunction(paraCoordinates(1), paraCoordinates(2));
					cartesianStressOnGivenPoint = shapeFuncs*cartesianStress;
					principalStress = ComputePrincipalStress(cartesianStressOnGivenPoint);					
					k1 = DirectionSelecting(iniDir, principalStress(typePSL), -principalStress(typePSL));	
					if 0 == AngleTerminationCondition(iniDir, k1), index = index-1; break; end
					%%k2
					midPot = nextPoint + k1*tracingStepWidth_/2;
					[elementIndex2, paraCoordinates2, bool1] = FindAdjacentElement(midPot);
					if ~bool1, index = index-1; break; end
					cartesianStress2 = cartesianStressField_(eNodMat_(elementIndex2,:)', :);
					shapeFuncs = ShapeFunction(paraCoordinates2(1), paraCoordinates2(2));
					cartesianStressOnGivenPoint2 = shapeFuncs*cartesianStress2;
					principalStress2 = ComputePrincipalStress(cartesianStressOnGivenPoint2);
					k2 = DirectionSelecting(k1, principalStress2(typePSL), -principalStress2(typePSL));		
					%%store	
					iniDir = k1;
					phyCoordList(index,:) = nextPoint;
					eleIndexList(index,:) = elementIndex;
					principalStressList(index,:) = principalStress;	
					%%next point
					nextPoint = nextPoint + tracingStepWidth_*k2;		
					[elementIndex, paraCoordinates, bool1] = FindAdjacentElement(nextPoint);
				end		
			end		
	end
	phyCoordList = phyCoordList(1:index,:);
	eleIndexList = eleIndexList(1:index,:);
	principalStressList = principalStressList(1:index,:);				
end

function val = AngleTerminationCondition(dirct1, dirct2)
	angle = acos((dirct1*dirct2') / (norm(dirct1)*norm(dirct2)));
	if angle > pi/6 %% continunity control
		val = 0;
	else
		val = 1;
	end
end

function targetDirection = DirectionSelecting(originalVec, Vec1, Vec2)
	normOriVec = norm(originalVec); normVec1 = norm(Vec1); normVec2 = norm(Vec2);
	angle1 = acos(originalVec*Vec1');
	angle2 = acos(originalVec*Vec2');
	if angle1 < angle2
		targetDirection = Vec1;
	else
		targetDirection = Vec2;
	end
end

function [nextElementIndex, paraCoordinates, opt] = FindAdjacentElement(physicalCoordinates)
	global nelx_; 
	global nely_; 
	global eleSize_;
	global nodeCoords_; 
	global eNodMat_;
	global carEleMapForward_;
	global boundingBox_;
	Lbound = boundingBox_(1,:);
	nextElementIndex = 0; paraCoordinates = []; opt = 0;
	
	physicalCoordinates = physicalCoordinates - Lbound;
	if 0==physicalCoordinates(1)
		eleX = 1;				
	else
		eleX = ceil(physicalCoordinates(1)/eleSize_);
		if eleX<1 || eleX>nelx_, return; end
	end
	if 0==physicalCoordinates(2)
		eleY = 1;
	else
		eleY = ceil(physicalCoordinates(2)/eleSize_);
		if eleY<1 || eleY>nely_, return; end
	end				
	
	tarEle = nely_*(eleX-1)+(nely_-eleY+1);
	nextElementIndex = carEleMapForward_(tarEle);
	if nextElementIndex	
		opt = 1;
		relatedNodes = eNodMat_(nextElementIndex,:);
		relatedNodeCoords = nodeCoords_(relatedNodes',:)-Lbound;
		paraCoordinates = 2*(physicalCoordinates - relatedNodeCoords(1,:)) / eleSize_ - 1;
	end
end

function ShowPSLs()
	global numEles_;
	global eNodMat_;
	global nodeCoords_;
	global nodesOnBoundary_;
	global majorPSLpool_;
	global minorPSLpool_;
	
	miniLength2Bshown = 5;
	for jj=1:length(majorPSLpool_)
		if majorPSLpool_(jj).length>miniLength2Bshown
			plot(majorPSLpool_(jj).phyCoordList(:,1), majorPSLpool_(jj).phyCoordList(:,2), ...
				'-', 'color', [252 141 98]/255, 'LineWidth', 3); hold on;		
		end
	end	
	for jj=1:length(minorPSLpool_)
		if minorPSLpool_(jj).length>miniLength2Bshown
			plot(minorPSLpool_(jj).phyCoordList(:,1), minorPSLpool_(jj).phyCoordList(:,2), ...
				'-', 'color', [102 194 165]/255, 'LineWidth', 3); hold on;
		end
	end	

	%%Show silhouette
	edgeIndices = eNodMat_(:, [1 2 2 1  2 3 3 2  3 4 4 3  4 1 1 4])';
	edgeIndices = reshape(edgeIndices(:), 4, 4*numEles_);	
	tmp = zeros(size(nodeCoords_,1),1); tmp(nodesOnBoundary_) = 1;
	tmp = tmp(edgeIndices); tmp = sum(tmp,1);
	boundaryEleEdges = edgeIndices(:,find(4==tmp));
	xPatchs = nodeCoords_(:,1); xPatchs = xPatchs(boundaryEleEdges);
	yPatchs = nodeCoords_(:,2); yPatchs = yPatchs(boundaryEleEdges);		
	cPatchs = zeros(size(yPatchs));
	hd = patch(xPatchs, yPatchs, cPatchs); hold on;
	set(hd, 'FaceColor', 'None', 'EdgeColor', [0.5 0.5 0.5], 'LineWidth', 2);
	
	axis equal; axis tight; axis off;
end

function CreatePSLsViaTSV(resCtrl)
	global nelx_;
	global nely_;
	global eleSize_;
	global boundingBox_;
	global nodeCoords_;
	global nodesOnBoundary_;
	global carNodMapBack_;
	
	global minimumEpsilon_;
	global mergeTrigger_;
	global seedPointsHistory_;
	global seedPoints_;
	global seedPointsValence_;
	global majorCoordList_; 
	global minorCoordList_;
	global startCoord_;
	global relaxedFactor_;
	
	global majorPSLpool_;
	global minorPSLpool_;

	%%1. Create Seed Points
	minimumEpsilon_ = min(boundingBox_(2,:)-boundingBox_(1,:))/resCtrl;
	validNodes = zeros((nelx_+1)*(nely_+1),1);
	validNodes(carNodMapBack_) = (1:length(carNodMapBack_))';
	validNodes = reshape(validNodes, nely_+1, nelx_+1);
	
	seedDensCtrl = max(ceil(minimumEpsilon_/eleSize_/5), 2);
	sampledNodes = validNodes(seedDensCtrl+1:seedDensCtrl:nely_+1-seedDensCtrl, ...
		seedDensCtrl+1:seedDensCtrl:nelx_+1-seedDensCtrl);	
	sampledNodes = reshape(sampledNodes, numel(sampledNodes), 1);
	sampledNodes(0==sampledNodes) = [];
	sampledNodes = setdiff(sampledNodes, nodesOnBoundary_);
	seedPointsHistory_ = nodeCoords_(sampledNodes,:);

	%%2. Initialize PSL pool
	startCoord_ = sum(boundingBox_,1)/2;
	majorCoordList_ = [];
	minorCoordList_ = [];
	mergeTrigger_ = minimumEpsilon_;
	seedPoints_ = seedPointsHistory_;
	numSeedPoints = size(seedPoints_,1);	
    seedPointsValence_ = zeros(numSeedPoints, 2);

	PreprocessSeedPoints();

	%%3. PSL sampling via TSV
	its = 0;
	looper = sum(sum(seedPointsValence_));	
	while looper<2*numSeedPoints
		its = its + 1;
		valenceMetric = sum(seedPointsValence_,2);
		unFinishedSpps = find(valenceMetric<2);
		spp = unFinishedSpps(1);
		if 0==looper
			[~, spp] = min(vecnorm(startCoord_-seedPoints_,2,2));
		else
			tmp	= seedPointsValence_(unFinishedSpps,:);
			tmp = find(1==sum(tmp,2));
			if ~isempty(tmp), unFinishedSpps = unFinishedSpps(tmp); end					
			[~, tarPos] = min(vecnorm(startCoord_-seedPoints_(unFinishedSpps,:),2,2));					
			spp = unFinishedSpps(tarPos);		
		end
		valences = seedPointsValence_(spp,:);
		seed = seedPoints_(spp,:);	
		if 0==valences(1)
			seedPointsValence_(spp,1) = 1;
			majorPSL = Have1morePSL(seed, 'MAJOR');		
			if 0==majorPSL.length
				looper = sum(sum(seedPointsValence_)); 
				disp([' Iteration.: ' sprintf('%4i',its) ' Progress.: ' sprintf('%6i',looper) ...
					' Total.: ' sprintf('%6i',3*numSeedPoints)]);
				continue; 
			end			
			majorPSLpool_(end+1,1) = majorPSL;				
			majorCoordList_(end+1:end+majorPSL.length,:) = majorPSL.phyCoordList;
			sppsEmptyMajorValence = find(0==seedPointsValence_(:,1));
			if ~isempty(sppsEmptyMajorValence)
				[potentialDisListMajor, potentialPosListMajor] = GetDisListOfPointList2Curve(seedPoints_(...
						sppsEmptyMajorValence,:), majorPSL.phyCoordList);					
				potentialSolidSppsMajor = find(potentialDisListMajor<relaxedFactor_);
				if ~isempty(potentialSolidSppsMajor)
					spps2BeMerged = sppsEmptyMajorValence(potentialSolidSppsMajor);
					seedPoints_(spps2BeMerged,:) = potentialPosListMajor(potentialSolidSppsMajor,:);								
					seedPointsValence_(spps2BeMerged,1) = 1;			
					modifiedMinorValences = HighCurvatureModification(spps2BeMerged, 'MINOR');				
					seedPointsValence_(modifiedMinorValences,2) = 1;	
				end
			end		
		end
		if 0==valences(2)
			seedPointsValence_(spp,2) = 1;			
			minorPSL = Have1morePSL(seed, 'MINOR');		
			if 0==minorPSL.length
				looper = sum(sum(seedPointsValence_)); 
				disp([' Iteration.: ' sprintf('%4i',its) ' Progress.: ' sprintf('%6i',looper) ...
					' Total.: ' sprintf('%6i',2*numSeedPoints)]);
				continue; 
			end		
			minorPSLpool_(end+1,1) = minorPSL;
			minorCoordList_(end+1:end+minorPSL.length,:) = minorPSL.phyCoordList;			
			sppsEmptyMinorValence = find(0==seedPointsValence_(:,2));
			if ~isempty(sppsEmptyMinorValence)   
				[potentialDisListMinor, potentialPosListMinor] = GetDisListOfPointList2Curve(seedPoints_(...
						sppsEmptyMinorValence,:), minorPSL.phyCoordList);					
				potentialSolidSppsMinor = find(potentialDisListMinor<relaxedFactor_);
				if ~isempty(potentialSolidSppsMinor)
					spps2BeMerged = sppsEmptyMinorValence(potentialSolidSppsMinor);
					seedPoints_(spps2BeMerged,:) = potentialPosListMinor(potentialSolidSppsMinor,:);
					seedPointsValence_(spps2BeMerged,2) = 1;				
					modifiedMajorValences = HighCurvatureModification(spps2BeMerged, 'MAJOR');					
					seedPointsValence_(modifiedMajorValences,1) = 1;					
				end
			end	
		end
		looper = sum(sum(seedPointsValence_));
		disp([' Iteration.: ' sprintf('%4i',its) ' Progress.: ' sprintf('%6i',looper) ...
			' Total.: ' sprintf('%6i',2*numSeedPoints)]);					
	end
end

function PreprocessSeedPoints()
	global seedPoints_;
	global seedPointsValence_;
	global majorPSLpool_; 
    global minorPSLpool_; 
	global relaxedFactor_;
	numMajorPSLs = length(majorPSLpool_);
	for ii=1:numMajorPSLs
		majorPSL = majorPSLpool_(ii);
		if majorPSL.length>0					
			sppsEmptyMajorValence = find(0==seedPointsValence_(:,1));
            if ~isempty(sppsEmptyMajorValence)
				[potentialDisListMajor, potentialPosListMajor] = GetDisListOfPointList2Curve(...	
					seedPoints_(sppsEmptyMajorValence,:), majorPSL.phyCoordList);
				potentialSolidSppsMajor = find(potentialDisListMajor<=relaxedFactor_);
				if ~isempty(potentialSolidSppsMajor)
					spps2BeMerged = sppsEmptyMajorValence(potentialSolidSppsMajor);							
					seedPoints_(spps2BeMerged,:) = potentialPosListMajor(potentialSolidSppsMajor,:);
					seedPointsValence_(spps2BeMerged,1) = 1;											
					modifiedMinorValences = HighCurvatureModification(spps2BeMerged, 'MINOR');
					seedPointsValence_(modifiedMinorValences,2) = 1;							
				end
			end
		end
	end
	numMinorPSLs = length(minorPSLpool_);
	for ii=1:numMinorPSLs
		minorPSL = minorPSLpool_(ii);
		if minorPSL.length>0	
			sppsEmptyMinorValence = find(0==seedPointsValence_(:,2));
            if ~isempty(sppsEmptyMinorValence)
				[potentialDisListMinor, potentialPosListMinor] = GetDisListOfPointList2Curve(...	
					seedPoints_(sppsEmptyMinorValence,:), minorPSL.phyCoordList);
				potentialSolidSppsMinor = find(potentialDisListMinor<=relaxedFactor_);
				if ~isempty(potentialSolidSppsMinor)
					spps2BeMerged = sppsEmptyMinorValence(potentialSolidSppsMinor);
					seedPoints_(spps2BeMerged,:) = potentialPosListMinor(potentialSolidSppsMinor,:);
					seedPointsValence_(spps2BeMerged,2) = 1;
					modifiedMajorValences = HighCurvatureModification(spps2BeMerged, 'MAJOR');
					seedPointsValence_(modifiedMajorValences,1) = 1;						
				end
			end
		end
	end
end

function iPSL = Have1morePSL(seed, psDir)
	global tracingStepWidth_;
	global boundingBox_;
	global snappingOpt_;
	stopCond = ceil(1.5*norm(boundingBox_(2,:)-boundingBox_(1,:))/tracingStepWidth_);	
	switch psDir
		case 'MAJOR'
			[iPSL, ~] = GeneratePrincipalStressLines(seed, [], stopCond, psDir);
		case 'MINOR'
			[~, iPSL] = GeneratePrincipalStressLines(seed, [], stopCond, psDir);
	end	
	if snappingOpt_, iPSL = CroppingPSLifNeeded(iPSL, psDir); end	
end

function [potentialDisList, potentialPosList] = GetDisListOfPointList2Curve(pointList, curveLine)
	global mergeTrigger_;
	disT = (curveLine(:,1) - pointList(:,1)').^2;
	disT = disT + (curveLine(:,2) - pointList(:,2)').^2;
	disT = sqrt(disT);	
	[minVal, minValPos] = min(disT,[],1);
	potentialDisList = minVal';
	potentialDisList = potentialDisList/mergeTrigger_;
	potentialPosList = curveLine(minValPos,:);	
end

function modifiedValences = HighCurvatureModification(spps2BeMerged, psDir)
	global majorCoordList_; global minorCoordList_;
	global seedPoints_;
	global seedPointsValence_;
	global mergeTrigger_;
	global relaxedFactor_;

	coordList = [];
	switch psDir
		case 'MAJOR'
			if isempty(majorCoordList_), modifiedValences = []; return; end
			coordList = majorCoordList_;
			spps2BeMerged = spps2BeMerged(find(0==seedPointsValence_(spps2BeMerged,1)));	
		case 'MINOR'
			if isempty(minorCoordList_), modifiedValences = []; return; end
			coordList = minorCoordList_;
			spps2BeMerged = spps2BeMerged(find(0==seedPointsValence_(spps2BeMerged,2)));
	end
	pointList = seedPoints_(spps2BeMerged,:);
	disT = (coordList(:,1) - pointList(:,1)').^2;
	disT = disT + (coordList(:,2) - pointList(:,2)').^2;
	disT = sqrt(disT);		
	minVal = min(disT, [], 1);
	minVal = minVal/mergeTrigger_;
	switch psDir
		case 'MAJOR'
			modifiedValences = find(minVal<relaxedFactor_);	
		case 'MINOR'
			modifiedValences = find(minVal<relaxedFactor_);	
	end	
	modifiedValences = spps2BeMerged(modifiedValences);
end


function tarPSL = CroppingPSLifNeeded(srcPSL, psDir)
	global mergeTrigger_;
	global majorCoordList_; 
    global minorCoordList_;
	tarPSL = srcPSL;
	if 5>=srcPSL.length, return; end
	disThreshold = 2;
	relaxedThreshold = 0.2;
	switch psDir
		case 'MAJOR'
			if isempty(majorCoordList_), return; end
			srcCoordList = majorCoordList_;						
		case 'MINOR'
			if isempty(minorCoordList_), return; end
			srcCoordList = minorCoordList_;		
	end
	
	if srcPSL.midPointPosition == srcPSL.length || srcPSL.midPointPosition == 1
		if 1==srcPSL.midPointPosition
			tarCoordList = tarPSL.phyCoordList;
			disT = (srcCoordList(:,1) - tarCoordList(:,1)').^2;
			disT = disT + (srcCoordList(:,2) - tarCoordList(:,2)').^2;	
			disT = sqrt(disT);
			miniDisList2SrcPSL = min(disT);
			tarPositions = find(miniDisList2SrcPSL<mergeTrigger_/disThreshold);
			if length(tarPositions)/size(tarCoordList,1)<relaxedThreshold, return; end
			if length(tarPositions) == size(tarCoordList,1), tarPSL = PrincipalStressLineStruct(); return; end
			startPos = 1; endPos = min(tarPositions);
		else
			tarCoordList = flip(tarPSL.phyCoordList,1);
			disT = (srcCoordList(:,1) - tarCoordList(:,1)').^2;
			disT = disT + (srcCoordList(:,2) - tarCoordList(:,2)').^2;	
			disT = sqrt(disT);
			miniDisList2SrcPSL = min(disT);
			tarPositions = find(miniDisList2SrcPSL<mergeTrigger_/disThreshold);
			if length(tarPositions)/size(tarCoordList,1)<relaxedThreshold, return; end
			if length(tarPositions) == size(tarCoordList,1), tarPSL = PrincipalStressLineStruct(); return; end
			startPos = srcPSL.length - min(tarPositions) + 1; endPos = srcPSL.length;
		end	
	else
		tarCoordList = tarPSL.phyCoordList(srcPSL.midPointPosition:srcPSL.length,:);
		disT = (srcCoordList(:,1) - tarCoordList(:,1)').^2;
		disT = disT + (srcCoordList(:,2) - tarCoordList(:,2)').^2;	
		disT = sqrt(disT);
		miniDisList2SrcPSL = min(disT);
		tarPositions = find(miniDisList2SrcPSL<mergeTrigger_/disThreshold);		
		if length(tarPositions)/size(tarCoordList,1)<relaxedThreshold
			endPos = srcPSL.length;
		elseif length(tarPositions) == size(tarCoordList,1)
			endPos = srcPSL.midPointPosition;
		else
			endPos = srcPSL.midPointPosition+min(tarPositions)-1;
		end
		
		tarCoordList = flip(tarPSL.phyCoordList(1:srcPSL.midPointPosition,:),1);
		disT = (srcCoordList(:,1) - tarCoordList(:,1)').^2;
		disT = disT + (srcCoordList(:,2) - tarCoordList(:,2)').^2;	
		disT = sqrt(disT);
		miniDisList2SrcPSL = min(disT);
		tarPositions = find(miniDisList2SrcPSL<mergeTrigger_/disThreshold);		
		if length(tarPositions)/size(tarCoordList,1)<relaxedThreshold
			startPos = 1;
		elseif length(tarPositions) == size(tarCoordList,1)
			startPos = srcPSL.midPointPosition;
		else		
			startPos = size(tarCoordList,1) - min(tarPositions) + 1;
		end
	end
	tarPSL.eleIndexList = srcPSL.eleIndexList(startPos:endPos,:);
	tarPSL.phyCoordList = srcPSL.phyCoordList(startPos:endPos,:);
	tarPSL.principalStressList = srcPSL.principalStressList(startPos:endPos,:);
	tarPSL.length = length(tarPSL.eleIndexList);
	tarPSL.midPointPosition = 1;		
end