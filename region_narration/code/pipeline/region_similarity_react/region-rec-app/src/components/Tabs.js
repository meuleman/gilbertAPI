import React from 'react';
import PropTypes from 'prop-types';
import Tab from './Tab';

class Tabs extends React.Component {
    static propTypes = {
        children: PropTypes.instanceOf(Array).isRequired,
    }

    constructor(props) {
        super(props);

        this.state = {
            activeTab: this.props.children[0].props.label,
        };
    }

    onClickedTabItem = (tab) => {
        this.setState({ activeTab: tab });
    }
}