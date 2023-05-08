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
    const chrVal = coordinates.split(':')[0].split('chr')[1]
    const start = coordinates.split(':')[1].split('-')[0]
    const end = coordinates.split(':')[1].split('-')[1]
    console.log()
    if (
      (coordinates) &&
      (chrVal.length > 0) &&
      (!isNaN(start)) &&
      (!isNaN(end)) && 
      (coordinates.slice(0, 3) === 'chr') && 
      (end > start) &&
      (this.state.chrValList.includes(chrVal))
    ) {
      this.props.setGenomicLocationProps(this.state.location)
      this.setState({
        location: ''
      })
    } else {
      alert('Please provide valid genomic coordinates')
    }
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
