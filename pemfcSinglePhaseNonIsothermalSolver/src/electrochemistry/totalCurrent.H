// calculating local current densities (Maricarmen Duarte)

// from
// E_cell = E_Nernst - I*ASR - (R*T/alphaAnode*F*ln(I*cMeRef/I0anode*concentrationMe)) 
//      - (R*T/alphaCathode*F*ln(I*cO2Ref/I0cathode*concentrationO2))
// E_cell = cell potential
// E_Nernst = Nernst potential
// I = local current density
// ASR = area specific resistance
// R = gas constant
// T = temperature
// alpha = transfer coefficient
// F = Faraday's constant
// I0anode = anode exchange current density
// I0cathode = cathode exchange current density


scalarField I(cathodePatch.size(), anodePatch.size());
{
    scalar A;
    scalar B;
    scalar r0 = upperBoundary.value();
    scalar r1 = lowerBoundary.value();

    forAll (I, facei)
    {
        A = (Rgas*cathodeT[facei])/(alphaCathodeElectrolyte.value()*F);
        B = (Rgas*anodeT[facei])/(alphaAnodeElectrolyte.value()*F);
 
               
        functionTotalCurrent ftc
        (
            Nernst[facei], 
            voltage.boundaryFieldRef()[electrolyteCathodeID][facei],
            R[facei],
            A,
            i0C[facei],
            B,
            i0A[facei],
            cO2Ref.value(),
            concO2[facei],
            cMeRef.value(),
            concMe[facei]
        );

        I[facei] = RiddersRoot<functionTotalCurrent>(ftc, 1.0e-5).root(r0, r1); // upper boundary for root bracketing must be less than the limiting current density

    }
}
