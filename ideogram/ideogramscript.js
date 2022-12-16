
  /* const fs = require("fs");
  const lines = fs.readFileSync("EpilogosScores.txt").toString().split("\r\n");

  let scores = [];

  for (let index = 0; index < lines.length; index+=1) {
   scores.push ({
      chromosome: lines[index]
    });
  }

  fs.writeFileSynch
  console.log("This is Epilogos Scores: " + lines);*/
  
  const CLPB = ["ISG15", "IFI6", "IFI44L", "IFI44", "PHGDH", "CMPK2", "RSAD2", "EIF2AK2", "MAL", "IFIH1", "SP100", "PARP9", "PLSCR1", "LAP3", "HERC6", "HERC5", "DDX60", "SAMD9L", "LY6E", "DDX58", "CHMP5", "ASS1", "IFIT2", "IFIT3", "IFIT1", "IFIT5", "IFITM2", "IFITM1", "IFITM3", "IRF7", "DRAP1", "CLPB", "LRRC32", "OAS1", "OAS3", "OAS2", "OASL", "EPSTI1", "IFI27", "CHAC1", "ISG20", "XAF1", "LGALS9", "DHX58", "IFI35", "PMAIP1", "BST2", "ZBP1", "HELZ2", "CBSL", "MX2", "MX1", "CBS", "USP18", "TYMP", "EGFL6", "CSAG3"];
  const RCL1 = ["ENO1", "MIR34AHG", "PIK3CD", "UBR4", "LDLRAP1", "AK4", "PHLDA3", "RRM2", "AAK1", "WDR54", "HSPD1", "PTMA", "NUP210", "XPC", "FAM162A", "RAB43", "ZMAT3", "ST6GAL1", "MFSD10", "SLC12A7", "CCNB1", "TCF7", "PTTG1", "HMGA1", "SESN1", "TRIM56", "GIMAP1", "MFHAS1", "BNIP3L", "RCL1", "STOML2", "ANXA1", "RGS3", "TNFSF8", "TMEM250", "TRIM22", "DDB2", "FTH1", "TIGAR", "TPI1", "TRIAP1", "CCDC92", "ITM2B", "ANXA2", "RPS27L", "MAP2K1", "PKM", "TMEM8A", "RSL1D1", "AC093512.2", "ALDOA", "ITGAL", "SMPD3", "KIF1C", "LASP1", "SMARCD2", "FDXR", "UBALD2", "TK1", "ACTG1", "KEAP1", "DNM2", "SMARCA4", "PRDX2", "PRKACA", "RPL18A", "LSM4", "SNRPD2", "FTL", "C19orf48", "PPP6R1", "SNPH", "SNRPB", "PI4KA", "MIF", "LGALS1", "SUN2", "APOBEC3H", "RPL3", "DENND6B", "PGK1", "MT-CO3"];
  const B3GALT2 = ["B3GALT2", "HACD2", "ACSS2"];
  const CHMP6 = ["RGS2", "DUSP10", "GUK1", "LYST", "WDR54", "IL18RAP", "CXCR4", "IFIH1", "KLHL24", "SQSTM1", "CD83", "TNFAIP3", "SNX9", "RALA", "INSIG1", "DUSP4", "ZBTB10", "SQLE", "DCTN3", "IFIT3", "NFKB2", "INPP5F", "SLC3A2", "NUDT22", "TPI1", "NELL2", "DDIT3", "DUSP6", "BTG1", "ISCU", "UBC", "NFKBIA", "PTGER2", "LGMN", "BMF", "PKM", "BCKDK", "AC004585.1", "GALK1", "CHMP6", "CYBC1", "TGIF1", "PMAIP1", "CNN2", "GPX4", "MYO1F", "ARRDC2", "LINC01480", "RELB", "PLA2G4C", "PPP1R15A", "IL4I1", "MIR663AHG", "ADA", "SREBF2", "TYMP", "GK", "CXCR3", "ABCD1"];
  const ZNF669 = ["MAP3K6", "ZNF669", "PKNOX1"];
  
  var annotsHeight = 3.5;
  fetchEnsembl = Ideogram.fetchEnsembl;

  var config = {
    organism: 'human',
    container: '#ideo-container',
    //chrHeight: 300,
  };


  function writeSelectedRange() {
    var r = ideogram.selectedRegion,
        from = r.from.toLocaleString(), // Adds thousands-separator
        to = r.to.toLocaleString(),
        extent = r.extent.toLocaleString();

    document.getElementById('from').innerHTML = from;
    document.getElementById('to').innerHTML = to;
    document.getElementById('extent').innerHTML = extent;
  }

  function writeCurrentPosition(position) {
    console.log("pos", position);
    document.getElementById('position').innerHTML = position.toLocaleString();
  }

  function displayAnnotsCaption() {
    console.log("enters display annots");
    var cfg = ideogram.config;
    var asm = cfg.assembly;
    var annotsUrl = cfg.annotationsPath;

    if (annotsUrl === '') {
      return;
    }

    d3.select('#annotations-summary').html(
      'Displaying ' + asm + ' BED file ' +
      '<a target="_blank" href="' + annotsUrl + '">' + annotsUrl + '</a>'
    )
    console.log("asm: " + asm + " annotsUrl: " + annotsUrl);
    console.log("rotate");
  }

  // This function executes ideogram after clicking on display BED
  d3.select('#display-bed-link').on('click', function() {
    var bedUrl = d3.select('#bed-url').nodes()[0].value;
    var asm = d3.select('input[name="asm"]:checked').nodes()[0].value;
    var anno = d3.select('input[name="anno"]:checked').nodes()[0].value;
    d3.select('#ideo-container').nodes()[0].innerHTML = '';
    ideogram = new Ideogram({
      organism: 'human',
      assembly: asm,
      annotationsPath: bedUrl,
      annotationHeight: annotsHeight,
      annotationsLayout: anno, // Switches from triangles to highlighted annotations
      resolution: 550,
      chrHeight: 300,
      chrWidth: 15,
      chrMargin: 2,
      orientation: 'vertical',
      rotatable: true,
      container: '#ideo-container',
      onDrawAnnots: displayAnnotsCaption
    });
  });


// START OF NEW IDEOGRAM TEST CODE


async function searchGenes(geneSymbols) {
  console.log("entered search genes");
  // Fetch positon of searched gene
  const organism = ideogram.config.organism;
  const orgUnderscored = organism.replace(/-/g, '_');
  let path = '/lookup/symbol/' + orgUnderscored;
  let body = {'symbols': geneSymbols}
  const ensemblGenes = await fetchEnsembl(path, body, 'POST');

  let annots = []

  ixnParam = interactions.map(ixn => `symbol:${ixn}`).join(' OR ')
  const taxid = ideogram.config.taxid;

  const myGeneBase = 'https://mygene.info/v3/query';
  const queryString = `?q=${ixnParam}&species=${taxid}&fields=symbol,genomic_pos`;
  let response = await fetch(myGeneBase + queryString);
  let data = await response.json();
  data.hits.forEach(hit => {
    annot = {
      name: hit.symbol,
      chr: hit.genomic_pos.chr,
      start: hit.genomic_pos.start,
      stop: hit.genomic_pos.end,
      id: hit.genomic_pos.ensemblgene,
      color: 'red'
    }
    annots.push(annot)
  })

  ideogram.drawAnnots(annots);

  document.querySelector('#ideo-container').style = 'display: none;';
  ideogram = new Ideogram({
    organism: 'human',
    assembly: 'GRCh38',
    annotations: [{
      h
    }],
    annotationHeight: annotsHeight,
    annotationsLayout: 'triangle',
    resolution: 550,
    chrHeight: 300,
    chrWidth: 15,
    chrMargin: 2,
    orientation: 'vertical',
    rotatable: true,
    container: '#ideo-container',
    onDrawAnnots: displayAnnotsCaption
  });
}

async function fetchInteractingGeneAnnots(interactions) {

  let annots = []

  ixnParam = interactions.map(ixn => `symbol:${ixn}`).join(' OR ')
  const taxid = ideogram.config.taxid;

  const myGeneBase = 'https://mygene.info/v3/query';
  const queryString = `?q=${ixnParam}&species=${taxid}&fields=symbol,genomic_pos`;
  let response = await fetch(myGeneBase + queryString);
  let data = await response.json();
  data.hits.forEach(hit => {
    annot = {
      name: hit.symbol,
      chr: hit.genomic_pos.chr,
      start: hit.genomic_pos.start,
      stop: hit.genomic_pos.end,
      id: hit.genomic_pos.ensemblgene,
      color: 'red'
    }
    annots.push(annot)
  })

  return annots;
}

async function fetchInteractingGenes(gene, organism) {
  let interactions = [];
  let orgNameSimple = organism.replace(/-/g, ' ');
  const queryString = `?query=${gene.name}&format=json`
  const url =
    `https://webservice.wikipathways.org/findInteractions${queryString}`;
  response = await fetch(url)
  data = await response.json()
  data.result.forEach(interaction => {
    if (interaction.species.toLowerCase() === orgNameSimple) {
      var right = interaction.fields.right.values;
      var left = interaction.fields.left.values;
      interactions = interactions.concat(right, left);
    }
  })

  // Filter out interactions that are definitely not gene symbols
  interactions = interactions.filter(ixn => {
    return ixn !== '' && !ixn.includes(' ') && ixn !== gene.name
  })

  interactions = Array.from(new Set(interactions)); // Deduplicate

  return interactions;
}



// END OF NEW IDEOGRAM TEST CODE

  // This function shows related genes on Ideogram
  function getGeneSelectedValue() {
    var geneSelected = document.getElementById('gene').value;
    if (geneSelected == "CLPB") {
      searchGenes(CLPB);

      /*ideogram = new Ideogram({
        organism: 'human',
        assembly: 'GRCh38',
        annotations: [{
          name: CLPB[1],
        }],
        annotationHeight: annotsHeight,
        annotationsLayout: 'triangle',
        resolution: 550,
        chrHeight: 300,
        chrWidth: 15,
        chrMargin: 2,
        orientation: 'vertical',
        rotatable: true,
        container: '#ideo-container',
        onDrawAnnots: displayAnnotsCaption
      });*/
    } else if (geneSelected == "RCL1") {
      searchGenes(RCL1);
      /*for (let i = 0; i < RCL1.length; i++) {
        searchGenes(RCL1[i]);
      }*/
      
    } else if (geneSelected == "B3GALT2") {
      searchGenes(B3GALT2);
    } else if (geneSelected == "CHMP6") {
      searchGenes(CHMP6);
    } else if (geneSelected == "ZNF669") {
      searchGenes(ZNF669);
    }

    console.log("Gene selected is: " + geneSelected);
  }


  // This function updates annotation styles
  d3.select('#update-annotations').on('click', function() {
    var bedUrl = d3.select('#bed-url').nodes()[0].value;
    var asm = d3.select('input[name="asm"]:checked').nodes()[0].value;
    var anno = d3.select('input[name="anno"]:checked').nodes()[0].value;
    var orient = d3.select('input[name="orient"]:checked').nodes()[0].value;
    d3.select('#ideo-container').nodes()[0].innerHTML = '';
    ideogram = new Ideogram({
      organism: 'human',
      assembly: asm,
      annotationsPath: bedUrl,
      annotationHeight: annotsHeight,
      annotationsLayout: anno, // Triangles to highlighted annotations
      resolution: 550,
      chrWidth: 15,
      chrMargin: 2,
      orientation: orient,
      rotatable: true,
      container: '#ideo-container',
      onDrawAnnots: displayAnnotsCaption
    });
  });

    // This function created more detailed annotations
    d3.select('#generate-detail').on('click', function() {
    var bedUrl = d3.select('#bed-url').nodes()[0].value;
    console.log("This is the bedURL: " + bedUrl)
    var asm = d3.select('input[name="asm"]:checked').nodes()[0].value;
    var anno = d3.select('input[name="anno"]:checked').nodes()[0].value;
    var orient = d3.select('input[name="orient"]:checked').nodes()[0].value;
    d3.select('#ideo-container-two').nodes()[0].innerHTML = '';

    function writeSelectedRange() {
      var r = ideogram.selectedRegion,
      from = r.from.toLocaleString(), // Adds thousands-separator
      to = r.to.toLocaleString(),
      extent = r.extent.toLocaleString();

      document.getElementById('from').innerHTML = from;
      document.getElementById('to').innerHTML = to;
      document.getElementById('extent').innerHTML = extent;
    }
    ideogram = new Ideogram({
      container: '.ideo-container-two',
      assembly: asm,
      organism: 'human',
      chromosome: '1',
      brush: 'chr1:104325484-119977655', // https://www.ncbi.nlm.nih.gov/dbvar/variants/nsv916356
      chrHeight: 900,
      resolution: 550,
      orientation: 'horizontal',
      onBrushMove: writeSelectedRange,
      onLoad: writeSelectedRange
    });
  });

  // This function generate a heatmap
  d3.select('#generate-heatmap').on('click', function() {
    var bedUrl = d3.select('#bed-url').nodes()[0].value;
    var asm = d3.select('input[name="asm"]:checked').nodes()[0].value;

    d3.select('#ideo-container').nodes()[0].innerHTML = '';
    var annotationTracks = [
      {id: 'expressionLevelTrack', displayName: 'Expression level'},
      {id: 'geneTypeTrack', displayName: 'Gene type'},
    ];

    var legend = [
      {
        name: 'Expression level',
        rows: [
          {color: '#88F', name: 'Low'},
          {color: '#CCC', name: 'Medium'},
          {color: '#F33', name: 'High'}
        ]
      },
      {
        name: 'Gene type',
        rows: [
          {color: '#00F', name: 'mRNA'},
          {color: '#0AF', name: 'misc_RNA'},
          {color: '#AAA', name: 'miRNA'},
          {color: '#FA0', name: 'tRNA'},
          {color: '#F00', name: 'lncRNA'}
        ]
      }
    ]

    var heatmaps = [
        {
          key: 'expression-level',
          thresholds: [
            ['2', '#88F'],
            ['4', '#CCC'],
            ['+', '#F33']]
        },
        {
          key: 'gene-type',
          thresholds: [
            ['0', '#00F'],
            ['1', '#0AF'],
            ['2', '#AAA'],
            ['3', '#FA0'],
            ['4', '#F00']
          ]
        }
      ]

    ideogram = new Ideogram ({
      container: '#ideo-container',
      organism: 'human',
      assembly: asm,
      resolution: 550,
      chrHeight: 275,
      chrWidth: 15,
      annotationHeight: annotsHeight,
      anootationsPath: bedUrl,
      annotationsLayout: 'heatmap',
      legend: legend,
      heatmaps: heatmaps,
      orient: 'vertical',
      annotationTracks: annotationTracks,
      rotatable: false, // Support for rotatable heatmaps is planned
      onDrawAnnots: displayAnnotsCaption
    });
  });

  // Handle uploaded BED file
  d3.select('#display-bed-upload').on('click', function() {
    let fileInput = document.getElementById('uploadInput');
    fileInput.onchange = () => {
        const bedObjects = [];
        const reader = new FileReader();
        
        reader.onload = (e) => {
          const data = e.target.result;
          tabDelimitedLines = data.split('\n').map((ln) => ln.split('\t'));
          tabDelimitedLines.forEach((elems) => {
            if (elems.length >= 3) {
              let bedObject = {
                chromosome: elems[0],
                start: parseInt(elems[1]),
                end: parseInt(elems[2]),
              };

              console.log("This is BED object: " + bedObject);
              if (elems.length >= 4) {
                bedObject = {...bedObject, id: elems[3]};
              }
              if (elems.length >= 5) {
                bedObject = {...bedObject, score: parseFloat(elems[4])};
              }
              if (elems.length >= 6) {
                bedObject = {...bedObject, strand: elems[5]};
              }
              bedObjects.push(bedObject);
            }
          });

            console.log(`${JSON.stringify(bedObjects, null, 2)}`);
            document.getElementById("uploadedInput").innerHTML=`<pre>${JSON.stringify(bedObjects, null, 2)}</pre>`;
        }
        
        for (let file of fileInput.files) {
            reader.readAsText(file);
        }
    }

    var asm = d3.select('input[name="asm"]:checked').nodes()[0].value;
    var anno = d3.select('input[name="anno"]:checked').nodes()[0].value;
    d3.select('#ideo-container').nodes()[0].innerHTML = '';
    ideogram = new Ideogram({
      organism: 'human',
      assembly: asm,
      annotations: [{
          name: 'test',
          chr: elems[0],
          start: parseInt(elems[1]),
          stop: parseInt(elems[2])
        }],
      annotationHeight: annotsHeight,
      annotationsLayout: anno, // Switches from triangles to highlighted annotations
      resolution: 550,
      chrHeight: 300,
      chrWidth: 15,
      chrMargin: 2,
      orientation: 'vertical',
      rotatable: true,
      container: '#ideo-container',
      onDrawAnnots: displayAnnotsCaption
    });

  });

  // Parse URL parameters
  var rawParams = document.location.search;
  console.log("RAWPARAMS = " + rawParams);
  var urlParams = {};
  var param, key, value;
  if (rawParams !== '') {
    rawParams = rawParams.split('?')[1].split('&');
    for (var i = 0; i < rawParams.length; i++) {
      param = rawParams[i].split('=');
      key = param[0];
      value = param[1];
      urlParams[key] = value;
    }
  }

  if ('annotsUrl' in urlParams) {
    // Example annotationsPath:
    // 'https://raw.githubusercontent.com/NCBI-Hackathons/Global_Screening_Arrays/master/files/output/pennTest2_gsrcCNVcall.bed',
    config.annotationsPath = urlParams.annotsUrl;
    config.annotationHeight = annotsHeight;
    config.onDrawAnnots = displayAnnotsCaption;
  }

  if ('asm' in urlParams) {
    config.assembly = urlParams.asm;
    d3.select('#' + urlParams.asm).nodes()[0].checked = true;
  }

  var ideogram = new Ideogram(config);