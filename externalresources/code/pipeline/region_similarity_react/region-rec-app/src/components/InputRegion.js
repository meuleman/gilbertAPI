import React from "react";

class InputRegion extends React.Component {
    state = {
        location: ""
    };

    onChange = (e) => {
        this.setState({
            [e.target.name]: e.target.value
        })
    };

    handleSubmit = (e) => {
        e.preventDefault();
        if (
            (this.state.location.trim()) && 
            (this.state.location.split(":")[0].split("chr")[1].length > 0) &&
            (!isNaN(this.state.location.split(":")[1].split("-")[0])) &&
            (!isNaN(this.state.location.split(":")[1].split("-")[1]))
        ) {
            this.props.setGenomicLocationProps(this.state.location);
            this.setState({
                location: ""
            })
        } else {
            alert("Please provide valid genomic coordinates")
        }
    };

    render() {
        return (
            <form onSubmit={this.handleSubmit} className="form-container">
                <input 
                    type="text"
                    className="input-text"
                    placeholder="Genomic Location"
                    value={this.state.location}
                    name="location"
                    onChange={(evt) => this.onChange(evt)}
                />
                <button className="input-submit">Go!</button>
            </form>
        )
    }
}

export default InputRegion