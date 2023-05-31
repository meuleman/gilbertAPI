function Spinner (props) {
  if (props.spinnerProps) {
  // if (true) {

    return (
      <div className="spinner-background">
        <div className='spinner'></div>
      </div>
      
    )
  } else {
    return (
      <div />
    )
  }
}
export default Spinner