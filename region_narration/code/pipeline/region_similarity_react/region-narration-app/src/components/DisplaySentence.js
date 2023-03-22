import React from 'react'
import factorsArr from '../factors.json'

function DisplaySentence (props) {
  if (props.sentenceProps.data !== '') {
    const queryRegion = props.sentenceProps.data[0]
    const targetFactors = queryRegion.target_factors
    const queryFactorOrder = queryRegion.query_factor_order
    const metrics = queryFactorOrder.map((ind) => {
      return factorsArr[ind].metric
    })
    const factors = queryFactorOrder.map((ind) => {
      return factorsArr[ind].name
    })
    const percentiles = queryRegion.percentiles
    const detailLevel = props.sentenceProps.detailLevel
    const thresh = Math.sqrt(0.9)
    let sentence = "This region is characterized by "
    const conversion = {
      'EnhA1': 'active enhancer 1',
      'EnhA2': 'active enhancer 2',
      'EnhBiv': 'bivalent enhancer',
      'EnhG1': 'genic enhancer 1',
      'EnhG2': 'genic enhancer 2',
      'EnhWk': 'weak enhancer',
      'Het': 'heterochromatin',
      'Quies': 'quiescent',
      'ReprPC': 'repressed polycomb',
      'ReprPCWk': 'weak repressed polycomb',
      'TssA': 'active TSS',
      'TssBiv': 'bivalent TSS',
      'TssFlnk': 'flanking TSS',
      'TssFlnkD': 'downstream flanking TSS',
      'TssFlnkU': 'upstream flanking TSS',
      'Tx': 'strong transcription',
      'TxWk': 'weak transcription',
      'ZNF/Rpts': 'ZNF genes + repeats',
      'Cancer/epithelial': 'cancer/epithelial',
      'Cardiac': 'cardiac',
      'Digestive': 'digestive',
      'Lymphoid': 'lymphoid',
      'Musculoskeletal': 'musculoskeletal',
      'Myeloid/erythroid': 'myeloid/erythroid',
      'Neural': 'neural',
      'Organdevel/renal': 'organ development/renal',
      'Placental/trophoblast': 'placental/trophoblast',
      'Primitive/embryonic': 'primitive/embryonic',
      'Pulmonarydevel': 'pulmonary development',
      'Renal/cancer': 'renal/cancer',
      'StromalA': 'stromal A',
      'StromalB': 'stromal B',
      'Tissueinvariant': 'tissue invariant',
      'Vascular/endothelial': 'vascular/endothelial'
    }

    let percentileCount = 0
    queryFactorOrder.map((regionIndex, index) => {
      if ((targetFactors.includes(regionIndex)) & (percentiles[index] > thresh)) {
        percentileCount += 1
      }
    })
    if (percentileCount == 0) {
      percentileCount = 1
    }
    const numFactorsToInclude = Math.min(percentileCount, detailLevel)

    let detailCount = 0
    queryFactorOrder.map((regionIndex, index) => {
      if ((targetFactors.includes(regionIndex)) & (detailCount < numFactorsToInclude)) {
        sentence += conversion[factors[index]] + " " + metrics[index]
        if ((numFactorsToInclude == 2) & (detailCount == 0)) {
          sentence += " and "
        } else if ((numFactorsToInclude > 2) & (detailCount == numFactorsToInclude - 2)) {
          sentence += ", and "
        } else if ((numFactorsToInclude > 2) & (detailCount < numFactorsToInclude - 2)) {
          sentence += ", "
        }
        detailCount += 1
      }
    })
    sentence += "."
    
    return (
      <div className='sentence-box'>
        <div className='sentence-header'><b>Region Summary</b></div>
        {/* <div className='sentence-body' dangerouslySetInnerHTML={{__html: props.sentenceProps.sentence}} /> */}
        <div className='sentence-body'>{sentence}</div>
      </div>
    )
  } else {
    return (
      <div />
    )
  }
}
export default DisplaySentence
