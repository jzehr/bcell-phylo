import React, { Component } from 'react';
import ReactDOM from 'react-dom';
import { BrowserRouter as Router, Route } from 'react-router-dom';
import { LinkContainer } from 'react-router-bootstrap';
import { Nav, Navbar, NavItem, NavDropdown, MenuItem, Grid } from 'react-bootstrap';
import 'bootstrap/dist/css/bootstrap.css';

import BCellPhylo from './bcell-phylo.js';
import JSONs from './jsons.js';

require('./main.css');



class App extends Component {
  constructor(props){
    super(props);
    this.patients = [ "28729", "48689", "67029", "77612", "78202", "93954", "99361", "99682", "GJS" ];
    this.genes = [
      "1-full", "1-18", "1-2", "1-24", "1-3", "1-46", "1-69", "1-69-2", "1-69D", "1-8", "1-45", "1-58",
      "2-full", "2-5", "2-70", "2-26",
      "3-full", "3-11", "3-13", "3-15", "3-20", "3-21", "3-23", "3-30-3", "3-33", "3-43", "3-48", "3-49", "3-53", "3-64", "3-66", "3-7", "3-72", "3-9", "3-43D", "3-73", "3-74",
      "4-full", "4-30-2", "4-30-4", "4-31", "4-34", "4-38-2", "4-39", "4-4", "4-59", "4-61",
      "5-full", "5-51", "5-10-1",
      "6-full", "6-1"
    ];
    this.state = { 
      patient: null,
      gene: null
    };
  }
  loadData(patient, gene, fragment) {
    this.setState({json: null}, function() {
      if (fragment == "full") {
        var json_path = `/data/${patient}/V${gene}.json`;
        console.log('full');
      } else {
        var json_path = `/data/${patient}/V${gene}-${fragment}.json`;
        console.log('fragment');
      }
      console.log(json_path);
      d3.json(json_path, (err, json_data) => {
        json_data.patient = patient;
        json_data.gene = gene;
        json_data.fragment = fragment;
        this.setState({
          patient: patient,
          gene: gene,
          fragment: fragment,
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
    const fragment_choices = {
      1: 18,
      2: 5,
      3: 11,
      4: '30-2',
      5: 51,
      6: 1
    };
    const patient = key.type == 'patient' ? key.value : this.state.patient;
    const gene = key.type == 'gene' ? key.value : this.state.gene;
    const fragment = key.type == 'fragment' ? key.value : fragment_choices[gene];
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
              {this.genes.filter(g=>g[0]==this.state.gene)
                .map(fragment => {
                  fragment = fragment.split('-').slice(1).join('-');
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
