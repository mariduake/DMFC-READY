    //cbpMesh: cathode side bipolarPlate mesh (J.P.Kone, 13.06.2016)

    fvMesh cbpMesh
    (
        IOobject
        (
            "cbp",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ
        )
    );

    // Read maps to global mesh
    // No mapping of face fields in solids.  HJ, 25/May/2011

    labelIOList cbpCellMap
    (
        IOobject
        (
            "cellRegionAddressing",
			cbpMesh.time().findInstance
            (
            	cbpMesh.meshDir(),
                "faces"
            ),
            polyMesh::meshSubDir,
			cbpMesh,
            IOobject::MUST_READ
        )
    );

    labelIOList cbpPatchesMap
    (
        IOobject
        (
            "boundaryRegionAddressing",
			cbpMesh.time().findInstance
            (
            	cbpMesh.meshDir(),
                "faces"
            ),
            polyMesh::meshSubDir,
			cbpMesh,
            IOobject::MUST_READ
        )
    );

