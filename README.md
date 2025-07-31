[IN DANISH]

Opsummering af scripts mappe i alfabetisk rækkefølge:

Autoplots: Genererer plots over de tz-gennemsnitslige profiler for dataen der behandles.

axialDistances: Script der fungerer som en R(s) funktion, og tillader at få den radialle position på et given sted på diverteren, relativt til strike-line.

create_data_folder & create_processed_folder: Benyttes til at lave data- og processed folders for nye simulationer, og håndterer hvis dataen allerede eksisterer etc., afhængig af hvilke options du giver dem.

dataIO: Simpelt objekt til at holde på data i et dictionary format, og benytter dumpProcessing til sin loading (alt gør)

densityShows: Nogle gamle ikke særlig kønne, men stadig informative plots der kan genereres af processed data, der viser densiteten af plasmaen over tid.

downstreamMap: Script der fungerer som en s(x) funktion, og tillader at mappe en outer midplane position i JT-60SA til dets korresponderende diverter position, ved at matche fluxen.

dumpProcessing: Kernen af alt datahåndtering. Dette er det første modul skrevet under projektet, og benyttes i mange sammenhænge. Adskillige scripts, herunder dataIO, heselProcessing, P_loss etc. benytter alle dumpProcessing direkte eller indirekte.

extraction: Gammelt script til data extraction fra før heselProcessing. En naiv (men simpel) implementering af at loade dataen. Er legacy kode, men har den stadig liggende i tilfælde af at den kunne blive brugbar. En ækvivalent funktionalitet kan findes i heselProcessing’s plainVariables method.

extraPlots: Indeholder plots for q (heat flux) med power fall-off length (PFOL), såvel som gamma (density flux) med density fall-off length (DFOL). Plots for energy over LCFS over tid er også inkluderet.

FOLs: Indeholder scripts til at beregne PFOL og DFOL.

fxEstimates: Script til at beregne et estimat for flux expansionen for et givet datasæt. Benytter PFOL.

heselProcessing: Dette script er en enkel class der benytter dataIO (og derfra dumpProcessing) til at loade alle de variable man kunne være interesseret i, i små bider, og er mit bedste forslag på en balancegang mellem CPU brug, og MEM brug. Hvis du skal implementere en ny variabel: Opdater listen over dependencies i starten af scriptet, tilføj beregningen som en method, og derefter tilføj den til calcVar metoden (den der henter dataen krævet for at beregne noget og derefter kalder metoden til at beregne størrelserne. Den ene kendte bug i øjeblikket, er at koden ikke kan håndtere hvis du giver den det samme x-index. Så hvis man f.eks. skal have noget omkring LCFS, så benytes LCFS og LCFS+1. Hvis du skal hente en variabel fra datafiler kan du benytte denne kode, specifikt plianVariables, averagedVariables (eller et af dets interfaces med t-, tz-, z- præfiks), eller summedVariables.

inputHandler: Håndterer inputfiler. Givet en dictionary med options henter den først inputfil-templaten fra template-mappen, og skriver en ny inputfil til en datamappe. Hvis der er options man ikke udfylder benytter den default options, og oplyser en om det. Simpel, og brugbar.

maxMeans: Inkluderer metode for at bestemme max energy, og 'max gennemsnitslig' energi, såvel som kode til at lave et plot over hvordan energien blev beregnet for et givet datasæt.

multiplot: Dette script genererer nogle af de samme plots som autoplots genererer, men i et samlet plot. Nice for at denne et hurtigt overblik over et run.

P_loss: Benytter heselProcessing til at ekstrahere q for elektronerne og ionerne, for at kunne beregne den samlede udadgående energi over LCFS. Har tre metoder, LCFS, maxMeans, og max. LCFS giver en konsekvent for lav værdi, så benyt en af de andre. Se mit speciale for detaljerede forklaringer.

processFolder: Meta-niveau funktion der indeholder en enkel funktion, som kalder heselProccesing på en datamappe med en masse variationer, man selv vælger. Ret nyttig, og benyttes i den samlede simulationsflow. Importeres i processing template.

q_projection: Script benyttet til at beregne q_D(s) - power profilen på diverten - for et givet datasæt. Fuld forklaring kan findes i speciale.

submissions: Core-script brugt til at udføre job submissions til Sophia clusteren. Hvis overført til en cluster der benytter et andet system end Sophia (slurm), skal dette opdateres. Samler alle de ovenstående scripts i en helhed med simsubmissions klassen, der inderholder do_full_process(), som kan tage en inputfil, sætte en simulation over med den, sætte et processing job over med dependency på simulationen, såvel som et visualiseringsjob med success dependency på processing jobbet. Ekstremt brugbar meta-level funktion. Bemærk at hvis man bare vil køre trin 2 og 3 og ikke 1 (sims), så er dette en mulighed.

submission_hack: Hurtigt skrevet script brugt til at køre submissions som i submissions filen, men hvor alt bliver submitted som et job i stedet for flere. Blev benyttet til at regenerere figurer uden at gøre folk et chock over mængden af jobs, og reducere overhead.

sweeps: Det primære interface til alt nHESEL processing. Tager et dictionary med options og laver et sweep over alle options givet, med rekursion for at generere alle mulige sammensætninger af de givne parametre (hvis ønsket), eller bare benytte default options for alle borset fra et parameter ad gangen (singulært sweep). Indeholder et par hjælpefunktioner for at gøre det lettere at håndtere sweeps.

sweeps_hack: Partnering til submissions_hack. Hurtigt skrevet script for at regenerere en masse figurere.
