/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Application
    testrhoCentralFoam2

Description
    Density-based compressible flow solver based on central-upwind schemes of
    Kurganov and Tadmor
    Generated from rhoCentralFoam of foamextend-3.0.1

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicFvMesh.H"   //new
#include "basicPsiThermo.H"
#include "zeroGradientFvPatchFields.H"
#include "fixedRhoFvPatchScalarField.H"
#include "motionSolver.H"    //new need check

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
#   include "setRootCase.H"

#   include "createTime.H"
//#   include "createMesh.H"     //this conflicts with createDynamicFvMesh
#   include "createDynamicFvMesh.H"   //new
#   include "createFields.H"
#   include "readThermophysicalProperties.H"
#   include "readTimeControls.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#   include "readFluxScheme.H"

    dimensionedScalar v_zero("v_zero",dimVolume/dimTime, 0.0);

    Info<< "\nStarting time loop\n" << endl;
    
    scalar CoNum = 0.0;
    scalar meanCoNum = 0.0;
    
    while (runTime.run())
    {
        #include "readTimeControls.H"
        #include "setDeltaT.H"
        
        runTime++;    // new position to update time

        Info<< "Time = " << runTime.timeName() << nl << endl;
        // --- upwind interpolation of primitive fields on faces
        // Do any mesh changes
       bool meshChanged =  mesh.update();


        surfaceScalarField rho_pos =
            fvc::interpolate(rho, pos, "reconstruct(rho)");
        surfaceScalarField rho_neg =
            fvc::interpolate(rho, neg, "reconstruct(rho)");

        surfaceVectorField rhoU_pos =
            fvc::interpolate(rhoU, pos, "reconstruct(U)");
        surfaceVectorField rhoU_neg =
            fvc::interpolate(rhoU, neg, "reconstruct(U)");

        volScalarField rPsi = 1.0/psi;
        surfaceScalarField rPsi_pos =
            fvc::interpolate(rPsi, pos, "reconstruct(T)");
        surfaceScalarField rPsi_neg =
            fvc::interpolate(rPsi, neg, "reconstruct(T)");

        surfaceScalarField e_pos =
            fvc::interpolate(e, pos, "reconstruct(T)");
        surfaceScalarField e_neg =
            fvc::interpolate(e, neg, "reconstruct(T)");

        surfaceVectorField U_pos = rhoU_pos/rho_pos;
        surfaceVectorField U_neg = rhoU_neg/rho_neg;

        surfaceScalarField p_pos = rho_pos*rPsi_pos;
        surfaceScalarField p_neg = rho_neg*rPsi_neg;

        surfaceScalarField phiv_pos = U_pos & mesh.Sf();
        surfaceScalarField phiv_neg = U_neg & mesh.Sf();
        Info<< "check point 1\n" << endl;
        
        Info << "last time step T: " << T.oldTime() << endl;
        Info << "last time step rho: " << rho.oldTime() << endl;

       // #include "volContinuity.H"   //newly added xiaoyue 
       /* 
        if (checkMeshCourantNo)
        {
            #include "meshCourantNo.H"
        }
*/
        
     //   Info << "mesh.phi() = " << mesh.phi() << endl;

        if (meshChanged)
	{
	  phiv_pos -= mesh.phi();
          phiv_neg -= mesh.phi();
          Info << "check mesh moving update1." << endl;
        }
       
// I think I need to check wether the code above has something to do with the temperature increasing


        //Info << "phiv_pos = " << phiv_pos << endl;
        //Info << "phiv_neg = " << phiv_neg << endl;
//floating point exception(core dumped)
/*       Info<< "check point 2\n" << endl;
       
       Info << "Cp:***************"<< endl
            << thermo.Cp() << endl;
       Info << "Cv:***************" << endl
            << thermo.Cv() << endl;
       Info << "rPsi:*************"<< endl  
	    << rPsi << endl;

       Info << "p:*************"<< endl  
	    << p << endl;
        Info << "rho:*************"<< endl  
	    << rho << endl;
*/
        Info << "cp 1.1" << endl;
        volScalarField c = sqrt(thermo.Cp()/thermo.Cv()*rPsi);
       
 
        Info << "cp 1.2" << endl;	
      
       Info<< "check point 2.1\n" << endl; 
   
        surfaceScalarField cSf_pos = fvc::interpolate(c, pos, "reconstruct(T)")*mesh.magSf();
        surfaceScalarField cSf_neg = fvc::interpolate(c, neg, "reconstruct(T)")*mesh.magSf();
       
        surfaceScalarField ap = max(max(phiv_pos + cSf_pos, phiv_neg + cSf_neg), v_zero);
        surfaceScalarField am = min(min(phiv_pos - cSf_pos, phiv_neg - cSf_neg), v_zero);
       Info<< "check point 2.2\n" << endl;
        surfaceScalarField a_pos = ap/(ap - am);
       Info<< "check point 2.3\n" << endl;
        surfaceScalarField amaxSf("amaxSf", max(mag(am), mag(ap)));
       Info<< "check point 2.4\n" << endl;
        surfaceScalarField aSf = am*a_pos;
        Info<< "check point 3 \n" << endl;
        if (fluxScheme == "Tadmor")
        {
            aSf = -0.5*amaxSf;
            a_pos = 0.5;
        }
        Info<< "check point 4\n" << endl;
        surfaceScalarField a_neg = (1.0 - a_pos);

        phiv_pos *= a_pos;
        phiv_neg *= a_neg;

        surfaceScalarField aphiv_pos = phiv_pos - aSf;
        surfaceScalarField aphiv_neg = phiv_neg + aSf;
    Info<< "check point 5\n" << endl;
        // Reuse amaxSf for the maximum positive and negative fluxes
        // estimated by the central scheme
        amaxSf = max(mag(aphiv_pos), mag(aphiv_neg));

        #include "compressibleCourantNo.H"
       

//        runTime++;      

        surfaceScalarField phi("phi", aphiv_pos*rho_pos + aphiv_neg*rho_neg);

        Info << "last time step aphiv_pos: " << aphiv_pos.oldTime() << endl;
        Info << "last time step rho_pos: " << rho_pos.oldTime() << endl;
        Info << "last time step aphiv_neg: " << aphiv_neg.oldTime() << endl;
        Info << "last time step rho_neg: " << rho_neg.oldTime() << endl;
// output flux variables
/*        Info << "FLUXES" << endl;
        Info << "aphiv_pos" << aphiv_pos << endl;
        Info << "rho_pos" << rho_pos << endl;
        Info << "aphi_neg" << aphiv_neg << endl;
        Info << "rho_neg" << rho_neg << endl;
*/  
     //Before solve density eqn, phi was defined here.
 //       Info << "phi before solve density eqn:" << phi << endl;

        surfaceVectorField phiUp =
            (aphiv_pos*rhoU_pos + aphiv_neg*rhoU_neg)
          + (a_pos*p_pos + a_neg*p_neg)*mesh.Sf();

        surfaceScalarField phiEp =
            aphiv_pos*(rho_pos*(e_pos + 0.5*magSqr(U_pos)) + p_pos)
          + aphiv_neg*(rho_neg*(e_neg + 0.5*magSqr(U_neg)) + p_neg)
          + aSf*p_pos - aSf*p_neg;
       Info<< "check point 6\n" << endl;
	// Make flux for pressure-work absolute

       if (meshChanged)
        {
            phiEp += mesh.phi()*(a_pos*p_pos + a_neg*p_neg);
            Info << "check Mesh moving update2." <<endl; 
        }
      
        Info<< "check point 7\n" << endl;
        volTensorField tauMC("tauMC", mu*dev2(fvc::grad(U)().T()));

  /*      Info << "Before solve density equation rho = " << rho << endl;
        Info << "Before solve density equation phi = " << phi << endl;
      
        Info << "Before solve density equation  p = " << p << endl;
        Info << "Before solve density equation  T = " << T << endl;
   */    
// --- Solve density
        solve(fvm::ddt(rho) + fvc::div(phi));
    /*    Info << "output updated rho values" << endl;
        Info << "rho = " << rho << endl;
        
        Info << "After solve density equation  p = " << p << endl;
        Info << "After solve density equation  T = " << T << endl;
*/
        Info << "last time step phi: " << phi.oldTime() << endl;
        Info << "last time step rho: " << rho.oldTime() << endl;

        // --- Solve momentum
        solve(fvm::ddt(rhoU) + fvc::div(phiUp));
       Info<< "check point 8\n" << endl;
        U.dimensionedInternalField() =
            rhoU.dimensionedInternalField()
           /rho.dimensionedInternalField();
        U.correctBoundaryConditions();
        rhoU.boundaryField() = rho.boundaryField()*U.boundaryField();

        volScalarField rhoBydt(rho/runTime.deltaT());
       Info<< "check point 9\n" << endl;
        if (!inviscid)
        {
            solve
            (
                fvm::ddt(rho, U) - fvc::ddt(rho, U)
              - fvm::laplacian(mu, U)
              - fvc::div(tauMC)
            );
            rhoU = rho*U;
        }

 //       Info << "After solve momentum eqn rho = " << rho << endl;
        //Info << "rhoE = " << rhoE << endl;
        // --- Solve energy
        surfaceScalarField sigmaDotU =
        (
            (
                fvc::interpolate(mu)*mesh.magSf()*fvc::snGrad(U)
              + (mesh.Sf() & fvc::interpolate(tauMC))
            )
            & (a_pos*U_pos + a_neg*U_neg)
        );
       Info<< "check point 10\n" << endl;
        solve
        (
            fvm::ddt(rhoE)
          + fvc::div(phiEp)
          - fvc::div(sigmaDotU)
        );
        Info << "check point 11\n" << endl;
        //Info << "rhoE = " << rhoE << endl;
        //Info << "After Energy eqn rho = " << rho << endl;
      // Info << "magSqr(U) = " << magSqr(U) <<endl;
        e = rhoE/rho - 0.5*magSqr(U);
        Info << "check 11.0" << endl;
        e.correctBoundaryConditions();
        
        Info << "check 11.1" << endl;

        thermo.correct();
        
        Info << "check 11.2" << endl;
        
        rhoE.boundaryField() =
            rho.boundaryField()*
            (
                e.boundaryField() + 0.5*magSqr(U.boundaryField())
            );
//        Info << "After energy eqn T = " << T << endl;


        if (!inviscid)
        {
            volScalarField k("k", thermo.Cp()*mu/Pr);
            solve
            (
                fvm::ddt(rho, e) - fvc::ddt(rho, e)
              - fvm::laplacian(thermo.alpha(), e)
              + fvc::laplacian(thermo.alpha(), e)
              - fvc::laplacian(k, T)
            );
            thermo.correct();
            rhoE = rho*(e + 0.5*magSqr(U));
        }

        p.dimensionedInternalField() =
            rho.dimensionedInternalField()
           /psi.dimensionedInternalField();
        p.correctBoundaryConditions();
        rho.boundaryField() = psi.boundaryField()*p.boundaryField();
        Info<< "check point 12\n" << endl;
        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
