/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    rhoSimpleFoam

Description
    Steady-state solver for turbulent flow of compressible fluids.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "fluidThermo.H"
#include "dynamicMomentumTransportModel.H"
#include "fluidThermophysicalTransportModel.H"
#include "simpleControl.H"
#include "pressureReference.H"
#include "fvModels.H"
#include "fvConstraints.H"
//有关等值面的头文件
#include "sampledIsoSurface.H"
#include "vtkSurfaceWriter.H"


// ************************************************************************* //
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "postProcess.H"

    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    #include "createFields.H"
    #include "createFieldRefs.H"
    #include "initContinuityErrs.H"
    //读取等值线字典
    dictionary isoSurfaceDict = IOdictionary(IOobject(
        "isoSurfaceDict",
        mesh.time().system(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE));
    //
    turbulence->validate();

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (simple.loop(runTime))
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        fvModels.correct();

        // Pressure-velocity SIMPLE corrector
        #include "UEqn.H"
        #include "EEqn.H"
        #include "pEqn.H"

        turbulence->correct();
        thermophysicalTransport->correct();
//
	cellSf=mesh.magSf();
	cellCenter=mesh.C();
	sumV=Foam::sum(cellVolume);
//      cellCenterX=mesh.C().component(0);
//	cellCenterY=mesh.C().component(1);		
//	cellCenterZ=mesh.C().component(2);

 //修改等值线字典中的isoValue数值
 //       isoSurfaceDict.set("isoValue", 0);
        //实例化isoSurface对象
        sampledSurfaces::isoSurface isosurf = sampledSurfaces::isoSurface(
            "isoSurface",
            mesh,
            isoSurfaceDict);
        //提取等值面
        isosurf.sample(cellCenterZ);
// or      isosurf.sample(mesh.C().component(2));
	//等值面的面单元向量
        faceList faces = isosurf.faces();
        //等值面的顶点坐标
        pointField points = isosurf.points();
        //计算等值面的面积
        scalar area = 0;
        forAll(faces, faceI)
        {
        cellSfiso=mag(faces[faceI].area(points));
            area += mag(faces[faceI].area(points));
        }
        
        Info << "iso面积:" << area << endl;
        //将等值面输出为VTK文件
        vtkSurfaceWriter vtkWriter = vtkSurfaceWriter(IOstream::streamFormat::ASCII);
        vtkWriter.write("postProcess",
                        "isoContours",
                        points,
                        faces);

	//
        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
