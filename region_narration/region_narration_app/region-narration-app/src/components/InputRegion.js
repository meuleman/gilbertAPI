import React from 'react'
// import alert from 'react-bootstrap/Alert';

class InputRegion extends React.Component {
  state = {
    location: '',
    chrValList: Array.from(
      {length: 22}, (value, index) => (index + 1).toString()
    ).concat(['X', 'Y'])
  }

  onChange = (e) => {
    this.setState({
      [e.target.name]: e.target.value
    })
  }

  handleSubmit = (e) => {
    e.preventDefault()
    const coordinates = this.state.location.trim()
    const coordinateFormat = /^chr\w+:\d+-\d+$/

    let passed = false
    if (coordinateFormat.test(coordinates)) {
      const coordSplit = coordinates.split(':')
      const chrVal = coordSplit[0].split('chr')[1]
      const start = coordSplit[1].split('-')[0]
      const end = coordSplit[1].split('-')[1]

      if (
        (end > start) &&
        (this.state.chrValList.includes(chrVal))
      ) {
        passed = true
      }
    }
    if (passed) {
      this.props.setGenomicLocationProps(this.state.location)
      this.setState({
        location: ''
      })
    } else {
      alert('Please provide valid genomic coordinates')
    }
    // console.log(coordinateFormat.test(coordinates))
  }

  render () {
    return (
      <form onSubmit={this.handleSubmit} className='form-container'>
        <input
          type='text'
          className='input-text'
          placeholder='Genomic Location'
          value={this.state.location}
          name='location'
          onChange={(evt) => this.onChange(evt)}
        />
        <button className='input-submit'>Go!</button>
      </form>
    )
  }
}

export default InputRegion
