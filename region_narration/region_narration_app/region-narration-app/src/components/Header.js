import React from 'react'

const Header = () => {
  const handleClick = () => {
    if (window.location.href !== document.location.href.split("?")[0]) {
      window.location.href = document.location.href.split("?")[0];
    }
  }
  return (
    <header>
      <h1 className='genomic-narration-header' onClick={() => handleClick()}>Genomic Narration</h1>
    </header>
  )
}

export default Header
