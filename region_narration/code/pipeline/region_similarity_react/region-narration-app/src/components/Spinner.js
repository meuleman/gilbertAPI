function Spinner (props) {
  if (props.spinnerProps) {
//   if (true) {

    return (
      <div className='spinner'></div>
    )
  } else {
    return (
      <div />
    )
  }
}
export default Spinner