Info<< "Reading thermophysical properties\n" << endl;
autoPtr<fluidThermo> pThermo
(
    fluidThermo::New(mesh)
);
fluidThermo& thermo = pThermo();
thermo.validate(args.executable(), "h", "e");

volScalarField& p = thermo.p();

//add
// 找到inlet这个patch的ID
label patchID = mesh.boundaryMesh().findPatchID("inlet"); 
// 创建inlet的polyPatch类型的引用
const polyPatch& cPatch = mesh.boundaryMesh()[patchID]; 
// 创建面标量场magSf，代表的是网格单元某一个面的面积
const surfaceScalarField& magSf = mesh.magSf(); 
// method1:
// 定义标量，并初始化为0.0
scalar patchArea = 0.0; 
// 遍历inlet上的所有单元面，然后将单元面的面积加起来
forAll(cPatch,faceI) 
{ 
    patchArea += magSf.boundaryField()[patchID][faceI]; 
} 
// 输出inlet这个边界面的面积
Info << "Patch Area " << patchArea << endl; 
// method2:
// 直接使用gSum函数来求inlet的面积
scalar AREA = gSum(magSf.boundaryField()[patchID]);
Info << "AREA = " << nl << AREA << endl;


//
Info<< "Reading field cellVolume\n" << endl;
volScalarField cellVolume
(
    IOobject
    (
        "cellVolume",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero",dimVolume,0)
);
cellVolume.ref()=mesh.V();
//
Info<< "Reading field sumV\n" << endl;
volScalarField sumV
(
    IOobject
    (
        "sumV",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);
//
Info<< "Reading field cellSf\n" << endl;
surfaceScalarField cellSf
(
    IOobject
    (
        "cellSf",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);
//
//
Info<< "Reading field cellSfiso\n" << endl;
volScalarField cellSfiso
(
    IOobject
    (
        "cellSfiso",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);
//
Info<< "Reading field cellCenter\n" << endl;
volVectorField cellCenter
(
    IOobject
    (
        "cellCenter",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);
//

//
Info<< "Reading field cellCenterX\n" << endl;
volScalarField cellCenterX
(
    IOobject
    (
        "cellCenterX",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

//cellCenterX=mesh.C(Foam::primitiveMesh::findNearestCell).component(0);
cellCenterX=mesh.C().component(0);
//
//
Info<< "Reading field cellCenterY\n" << endl;
volScalarField cellCenterY
(
    IOobject
    (
        "cellCenterY",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);
cellCenterY=mesh.C().component(1);
//
//
Info<< "Reading field cellCenterZ\n" << endl;
volScalarField cellCenterZ
(
    IOobject
    (
        "cellCenterZ",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);
cellCenterZ=mesh.C().component(2);
//
volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    thermo.rho()
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

#include "compressibleCreatePhi.H"

pressureReference pressureReference
(
    p,
    simple.dict(),
    thermo.incompressible()
);

mesh.setFluxRequired(p.name());

Info<< "Creating turbulence model\n" << endl;
autoPtr<compressible::momentumTransportModel> turbulence
(
    compressible::momentumTransportModel::New
    (
        rho,
        U,
        phi,
        thermo
    )
);

Info<< "Creating thermophysical transport model\n" << endl;
autoPtr<fluidThermophysicalTransportModel> thermophysicalTransport
(
    fluidThermophysicalTransportModel::New(turbulence(), thermo)
);

dimensionedScalar initialMass = fvc::domainIntegrate(rho);

#include "createMRF.H"
#include "createFvModels.H"
#include "createFvConstraints.H"
