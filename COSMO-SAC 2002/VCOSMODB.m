function [Vcosmo] = VCOSMODB(comp)

for i = 1:size(comp,2)
    component = comp(i);
    switch component
        case "water"
            Vcosmo(i) = 25.73454;
        case "ethanol"
            Vcosmo(i) = 70.19948;
        case "butanol"
            Vcosmo(i) = 113.87358;
        case "methanol"
            Vcosmo(i) = 48.77104;
        case "propanol"
            Vcosmo(i) = 92.10607;
        case "2propanol"
            Vcosmo(i) = 92.41349;
        case "2butanol"
            Vcosmo(i) = 113.74155;
        case "isobutanol" % 2-methyl-1-propanol
            Vcosmo(i) = 113.93402;
        case "tertbutanol" % 2-methyl-2-propanol
            Vcosmo(i) = 114.85208;
        case "pentanol"
            Vcosmo(i) = 135.91141;
        case "2pentanol"
            Vcosmo(i) = 135.06275;
        case "aceticacid"
            Vcosmo(i) = 74.68827;
        case "ethylacetate"
            Vcosmo(i) = 118.1693;
        case "butylacetate"
            Vcosmo(i) = 162.13696;
        case "isobutylacetate"
            Vcosmo(i) = 162.83181;
        case "isoamylacetate"
            Vcosmo(i) = 184.38148;
        case "isopentylacetate"
            Vcosmo(i) = 184.38148;
        case "butyricacid"
            Vcosmo(i) = 117.70778;
        case "heptylacetate"
            Vcosmo(i) = 227.62932;
        case "butanal"
            Vcosmo(i) = 108.16459;
        case "diisopropylether"
            Vcosmo(i) = 159.52489;
        case "hexane"
            Vcosmo(i) = 146.12927;
        case "heptane"
            Vcosmo(i) = 167.98102;
        case "heptene"
            Vcosmo(i) = 163.39637;
        case "oxylene"
            Vcosmo(i) = 153.88669;
        case "mxylene"
            Vcosmo(i) = 153.76642;
        case "pxylene"
            Vcosmo(i) = 154.46557;
        case "mesitylene" % 1,3,5-trimethylbenzene
            Vcosmo(i) = 175.55351;
        case "benzene"
            Vcosmo(i) = 110.22176;
        case "14dioxane"
            Vcosmo(i) = 110.13815;
    end        
end