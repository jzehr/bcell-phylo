import React, { Component } from 'react';
import ReactDOM from 'react-dom';
import { BrowserRouter as Router, Route } from 'react-router-dom';
import { LinkContainer } from 'react-router-bootstrap';
import { Nav, Navbar, NavItem, NavDropdown, MenuItem, Grid } from 'react-bootstrap';
import 'bootstrap/dist/css/bootstrap.css';

import BCellPhylo from './bcell-phylo.js';
import JSONs from './jsons.js';

require('./main.css');
const patient_v_pairs = require('../data/patient_v_pairs.json'); 



class App extends Component {
  constructor(props){
    super(props);
    this.patients = Object.keys(patient_v_pairs);
    this.state = { 
      patient: null,
      fragments: [],
      fragment: null,
      genes: [],
      gene: null
    };
  }
  loadData(patient, gene, fragment) {
    this.setState({json: null}, function() {
      if (fragment == "full") {
        var json_path = `/data/${patient}/V${gene}.json`;
      } else {
        var json_path = `/data/${patient}/V${gene}-${fragment}/dashboard.json`;
      }
      d3.json(json_path, (err, json_data) => {
        json_data.patient = patient;
        json_data.gene = gene;
        json_data.fragment = fragment;
        this.setState({
          patient: patient,
          gene: gene,
          genes: Object.keys(patient_v_pairs[patient]),
          fragment: fragment,
          fragments: patient_v_pairs[patient][gene],
          json: json_data
        });
      });
    });
  }
  componentDidMount() {
    const patient = '77612';
    const gene = '3';
    const fragment = '11';
    this.loadData(patient, gene, fragment);
  }
  onSelect(key){
    var patient, gene, fragment;
    if(key.type == 'patient') {
      patient = key.value;
      gene = Object.keys(patient_v_pairs[patient])[0];
      fragment = patient_v_pairs[patient][gene][0];
    } else if(key.type == 'gene') {
      patient = this.state.patient;
      gene = key.value;
      fragment = patient_v_pairs[patient][gene][0];
    } else if(key.type == 'fragment') {
      patient = this.state.patient;
      gene = this.state.gene;
      fragment = key.value;
    }
    this.loadData(patient, gene, fragment);
  }
  render(){
    return(<Router>
      <div>
        <Navbar onSelect={key=>this.onSelect(key)} fixedTop>
          <Navbar.Header>
            <Navbar.Brand>
              ACME Immunology Visualization
            </Navbar.Brand>
          </Navbar.Header>
          <Nav>
            <NavDropdown title='Patient' id='patient'>
              {this.patients.map(patient_id => {
                const eventKey = { type: 'patient', value: patient_id };
                return (<MenuItem
                  eventKey={eventKey}
                  active={this.state.patient == patient_id}
                  key={patient_id}
                >
                  {patient_id}
                </MenuItem>);
                }) }
            </NavDropdown>
            <NavDropdown title='Gene' id='patient'>
              {[1, 2, 3, 4, 5, 6].map(gene => {
                const eventKey = { type: 'gene', value: gene };
                return (<MenuItem
                  eventKey={eventKey}
                  active={this.state.gene == gene}
                  key={gene}
                >
                  {gene}
                </MenuItem>);
              }) }
            </NavDropdown>
            <NavDropdown title='Fragment' id='fragment'>
              {this.state.fragments
                .map(fragment => {
                  const eventKey = { type: 'fragment', value: fragment };
                  return (<MenuItem
                    eventKey={eventKey}
                    active={this.state.fragment == fragment}
                    key={this.state.gene+'-'+fragment}
                  >
                    {fragment}
                  </MenuItem>);
                }) }
            </NavDropdown>
          </Nav>
        </Navbar>
        <Grid fluid>
          <Route exact path="/" render={ () => <BCellPhylo json={this.state.json} /> } />
        </Grid>
      </div>
    </Router>);
  }
}

ReactDOM.render(
  <App />,
  document.body.appendChild(document.createElement('div'))
)
