# PALEOocean Reactions

## Ocean geometry and circulation transport
```@meta
CurrentModule = PALEOocean.Ocean
```

```@docs
OceanNoTransport.ReactionOceanNoTransport
OceanTransport3box.ReactionOceanTransport3box
OceanTransport6box.ReactionOceanTransport6box
OceanTransportColumn.ReactionOceanTransportColumn
OceanTransportRomaniello2010.ReactionOceanTransportRomanielloICBM
OceanTransportTMM.ReactionOceanTransportTMM
```

## Vertical Transport
```@docs
VerticalTransport.ReactionLightColumn
VerticalTransport.ReactionExportDirect
VerticalTransport.ReactionExportDirectColumn
VerticalTransport.ReactionSinkFloat
```

## Biological Production
```@docs
BioProd.ReactionBioProdPrest
BioProd.ReactionBioProdMMPop
```

## Ocean surface air-sea flux
```@meta
CurrentModule = PALEOocean.Oceansurface
```

```@docs
AirSeaExchange.ReactionAirSea
AirSeaExchange.ReactionAirSeaO2
AirSeaExchange.ReactionAirSeaCO2
AirSeaExchange.ReactionAirSeaCH4
AirSeaExchange.ReactionAirSeaFixedSolubility
```

## Ocean floor burial

```@meta
CurrentModule = PALEOocean.Oceanfloor
```


###  Carbonate burial
```@docs
Burial.ReactionShelfCarb
Burial.ReactionBurialEffCarb
```

### Organic carbon and phosphorus burial
```@docs
Burial.ReactionBurialEffCorgP
```

## Global
```@meta
CurrentModule = PALEOocean.Global
```
```@docs
Insolation.ReactionForceInsolationModernEarth
Insolation.insolMITgcmDIC
```
