/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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
    pisoFoam

Description
    Transient solver for incompressible, turbulent flow, using the PISO
    algorithm.

    Sub-models include:
    - turbulence modelling, i.e. laminar, RAS or LES
    - run-time selectable MRF and finite volume options, e.g. explicit porosity

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "pisoControl.H"
#include "fvOptions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "postProcess.H"

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
	#include "createTimeControls.H"
    #include "createControl.H"
    #include "createFields.H"
    #include "createFvOptions.H"
    #include "initContinuityErrs.H"

    turbulence->validate();

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.loop())
    {

		//---Turn CT# data into porosity and permeability data, but only during first timestep after 0
		if ( runTime.value() == runTime.deltaT().value() )	
			{
				forAll (mesh.cells(),celli)
				{
				if ( CT[celli] >= CTupper.value() )
					 { eps[celli] = epsLower.value();}
				else 
					if ( CT[celli] <= CTlower.value() )
					{  eps[celli] = epsUpper.value();}
					else { eps[celli] = (((epsLower.value() - epsUpper.value())
									  * (CT[celli] - CTlower.value()))
									  / (CTupper.value() - CTlower.value())) + epsUpper.value(); }

//				K[celli] =	K0.value() * Foam::pow((1-eps[celli]),2.0) / Foam::pow(eps[celli],3.0);
			   	K[celli] = K0.value()
						 * (Foam::pow(1-eps[celli],2.0)
						 / Foam::pow(max(epsLower.value(),eps[celli]-epsCrit.value()),Kn.value()));

				}
			Info<< "eps & K set" << nl << endl;
			}

		 Info<< "Time = " << runTime.timeName() << nl << endl;

        #include "readTimeControls.H"
        #include "CourantNo.H"
		#include "setDeltaT.H"
			

        //---Pressure-velocity PISO corrector
        {
            #include "UEqn.H"

            // --- PISO loop
            while (piso.correct())
            {
                #include "pEqn.H"
            }
        
		}

        laminarTransport.correct();
        turbulence->correct();
		Uint = U/eps;

		//Limit bacteria attachment and CaCO3 precipitation to surfaces
		const scalar poreScale(readScalar(transportProperties.lookup("poreScale")));
		if( poreScale == 1)
		{
			//---Define attachment region (6-connected). N.B. not across processor boundaries
		    forAll(mesh.cells(),cellI)
		    {
				//define solids as any cell with porosity below threshold
				if (eps[cellI] <= solidThreshold.value())
				{
					solid[cellI] = 1.0;
				}	
				//for each cell, attach value is sum of neighbouring cells solid value
				labelList adjacent = mesh.cellCells()[cellI];
				int size = adjacent.size();
				for(int j=0; j<size; j++)
				{			
					attach[cellI] += solid[adjacent[j]];			
				}
				//substitute sum of solid value with 1, indicating attachment can occur
				//in any cell with at least one neighbouring solid cell
				if (attach[cellI] > 0)
				{
					attach[cellI] = 1.0;
				}
				//Additional attachment seeding independent of porosity, specified through CT# field
				if (CT[cellI] >= attachSeedThreshold.value())
				{
					attach[cellI] = 1.0;
				}
			}
		}
		else 
		    forAll(mesh.cells(),cellI)
		    {
				attach[cellI] = 1.0;
			}


		// --- Determine velocity related bacteria attachment --- //
		forAll (mesh.cells(),celli)
		{
		if ( mag(Uint[celli]) >= zAttachUpperThreshold.value() )
			 { zAttachSwitch[celli] = scalar(0.0); }
		else 
			if ( mag(Uint[celli]) <= zAttachLowerThreshold.value() )
			 { zAttachSwitch[celli] = scalar(1.0); }
			 else
			 { zAttachSwitch[celli] = scalar(1.0) - ((mag(Uint[celli]) 
									- zAttachLowerThreshold.value()) 
									/ (zAttachUpperThreshold.value() 
									- zAttachLowerThreshold.value())) ; }
		}

		// --- bacteria transport, 		    
        forAll(mesh.cells(),cellI)
        {
            //define solids as any cell with porosity below threshold
            if (eps[cellI] <= solidThreshold.value())
            {
                zGrowth[cellI] = 0.0;
            }	
            if (eps[cellI] > solidThreshold.value())
            {
                zGrowth[cellI] = 1.0;
            }
        }
        solve
        (
            fvm::ddt(eps, bacteria)
          + fvm::div(phi, bacteria)
          - fvm::laplacian((eps*D_bacteria), bacteria)
		  + fvm::Sp(((zAttach * zAttachSwitch) + zStrain) * attach * eps * zGrowth, bacteria)
			// - fvm::Sp(YCOD * (zBacteria+bacteria) * qhat_COD * eps * (PO4/(KPO4+PO4)) * (NO3/(KNO3+NO3)) * (1/(KCOD+COD)), COD)
			- fvm::Sp(YCOD * qhat_COD * eps * (PO4/(KPO4+PO4)) * (NO3/(KNO3+NO3)) * (COD/(KCOD+COD)), bacteria)
        );

		// --- bacteria attachment, decay, encapsulation --- //
        solve
        (
            fvm::ddt(eps, zBacteria)
      	  - fvc::Sp(((zAttach * zAttachSwitch) + zStrain) * attach * eps * zGrowth, bacteria)
      	  + fvm::Sp(zDecay, zBacteria)
          + fvm::Sp(YCOD * qhat_COD * eps * (PO4/(KPO4+PO4)) * (NO3/(KNO3+NO3)) * (COD/(KCOD+COD)) * zGrowth, zBacteria)
     	  // + fvm::Sp(fvc::ddt(CaCO3) / zEncapsulationConst, zBacteria) //should porosity be in here?
        );
        // aaa = 
        // Info<< "test = " << min(mag(YCOD * qhat_COD * eps * (PO4/(KPO4+PO4)) * (NO3/(KNO3+NO3)) * (COD/(KCOD+COD)) * (1-eps/0.2))) << " s" << endl;
		// --- COD transport, and consumption --- //
        solve
        (
            fvm::ddt(eps, COD)
          + fvm::div(phi, COD)
          - fvm::laplacian((eps*D_COD), COD)
		  + fvm::Sp((zBacteria+bacteria) * qhat_COD * eps * (PO4/(KPO4+PO4)) * (NO3/(KNO3+NO3)) * (1/(KCOD+COD)), COD)
        );

		// --- nitrate transport and ureolysis source --- //
        solve
        (
            fvm::ddt(eps, NO3)
          + fvm::div(phi, NO3)
          - fvm::laplacian((eps*D_NO3), NO3)
          + fvm::Sp((YCOD/YNO3) * (zBacteria+bacteria) * qhat_COD * eps * (PO4/(KPO4+PO4)) * (COD/(KCOD+COD)) * (1/(KNO3+NO3)), NO3)
        );

		// --- phosphate transport, consumption --- //
        solve
        (
            fvm::ddt(eps, PO4)
          + fvm::div(phi, PO4)
          - fvm::laplacian((eps*D_PO4), PO4)
		  + fvm::Sp((YCOD/YPO4) * (zBacteria+bacteria) * qhat_COD * eps * (COD/(KCOD+COD)) * (NO3/(KNO3+NO3)) * (1/(KPO4+PO4)), PO4)
        );

		// // -- explicit Michaelis-Menten hydrolysis --- //
		// 	hydrolysis = fvc::Sp((zBacteria+bacteria) * Ku * eps * (1/(Km+COD)) * (KNO3/(KNO3+NO3)), COD);

		// // --- calcium transport and precipitation sink --- //
	  //   solve
    //     (
    //         fvm::ddt(eps, Ca)
    //       + fvm::div(phi, Ca)
    //       - fvm::laplacian((eps*D_Ca), Ca)
		//   + fvc::Sp(Kp * min(Ca,PO4) * eps, attach)
    //     );


	
		// --- porosity and permeabiltiy change due to attached bacteria --- //
        solve
        (
            fvm::ddt(eps)
          + fvc::ddt(zBacteria) / rhozBacteria
        );
        // Info<< "test = " << min(eps) << " s"
        //     << "  ClockTime = " << max(mag(U)) << "," <<  min(mag(U)) << " aa" << max(mag(phi)) << " s"
        //     << nl << endl;
		K =	K0 * pow((1-eps),2.0) / pow(max(epsLower.value(),eps-epsCrit),Kn);

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
