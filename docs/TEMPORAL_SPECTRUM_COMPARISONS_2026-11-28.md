# Four saved-spectrum comparisons and nine-frame geometry

28 November 2026. [Full visual report](../results/temporal_spectrum_comparisons/20261128T141633Z_2ac2c0f8/report.md), [timing diagram](../results/temporal_spectrum_comparisons/20261128T141633Z_2ac2c0f8/window_geometry.png), [raw/Zoverview](../results/temporal_spectrum_comparisons/20261128T141633Z_2ac2c0f8/four_cases_spectral_evolution.png).

Input4096samples/4s;added tone(s)1024samples,support[1.25,2.25)s. W1Hann1024samples,hop128samples=.125s gives25completewindows,centers.5..3.5s. Central9are chosen by center1.25..2.25inclusive,j6..14;their input windows run[.75,1.75)through[1.75,2.75),union2s. Only center1.75covers the entire one-second nominal pulse. Fifteenfull-mapwindows overlap the pulse at all;central9are not all signal-touching windows. Overlap87.5%,notindependent measurements;Hannweights change the contributions.

Four fixed illustrative cases,allresultantSNR18: true88/90equal passes5/9,maxrun3;matchedratio.5passes2/9,maxrun2;wrong86/89fortrue86/88passes5centralconsecutiveframes;single102producesanapparent102/105pairin3centralconsecutiveframes. The full report shows9successive normalized andrawspectra percase,not just peak heights. Raw/normalizedspectralmapsuse sharedcolorrangesacrosscases;profileaxes fixedwithin eachcase. Allnumeric localspectralvalues andwindowgeometryexported.

Attrueequalpair’s center1.75,peaks30.09/29.02arehighbutvalley0.940fails;theweakertruepaircanremainaboveheight4whilehavingnoreliabledip. Wrong/singlecompanionscanlookstableandisolated. This suggests inspecting fullshape/jointevolution,butdoesnotestablisharule. Firsttwocasesarepaired;othersaredifferentbases. Selectionisillustrative,notheld-outvalidation. Insufficiencyofthecurrentstatisticsisnotproofthatnoothermethodcanusewaveforminformation.

Parent529andtemporal21artifacthashesunchanged;36gateand25geometrychecksPASS. No waveform/mapgeneration,detector/gatechanges,O3,commit/push. Newgalleryandfindingsarelocalforreview.
