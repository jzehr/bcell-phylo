import React, { Component } from 'react';
import ReactDOM from 'react-dom';
import { BrowserRouter as Router, Route, withRouter } from 'react-router-dom';
import { LinkContainer } from 'react-router-bootstrap';
import { Nav, Navbar, NavItem, NavDropdown, MenuItem, Grid } from 'react-bootstrap';
import 'bootstrap/dist/css/bootstrap.css';

import BCellInterface from './bcell-phylo.js';
import Overview from './overview.js';

require('./main.css');
const patient_v_pairs = require('../data/patient_v_pairs.json'); 
const patients = Object.keys(patient_v_pairs);


function parse_route(route) {
  const split_route = route.split(/&|-|\//),
    parsed_route = {
      patient_id: +split_route[1] || 77612,
      gene: split_route[2] ? String(split_route[2])[1] : "3",
      allele: split_route[3] ? String(split_route[3]) : "11"
    };
  return parsed_route;
}

function get_route(patient_id, gene, allele) {
  return `/${patient_id}&V${gene}-${allele}`;
}

const PatientDropdown = withRouter(function(props) {
  return (<NavDropdown title='Patient' id='patient'>
    {patients.map(patient_id => {
      const gene = 3,
        allele = patient_v_pairs[patient_id][String(gene)][0],
        route = get_route(patient_id, gene, allele);
      return (<LinkContainer to={route} key={patient_id}>
        <MenuItem>
          {patient_id}
        </MenuItem>
      </LinkContainer>);
    }) }
  </NavDropdown>);
});

const FamilyDropdown = withRouter(function(props) {
  const { location } = props,
    { patient_id } = parse_route(location.pathname);
  return (<NavDropdown title='Family' id='family'>
    {[1, 2, 3, 4, 5, 6].map(gene => {
      const allele = patient_v_pairs[patient_id][gene][0],
        route = get_route(patient_id || 77612, gene, allele);
      return (<LinkContainer to={route} key={gene}>
        <MenuItem>
          {gene}
        </MenuItem>
      </LinkContainer>);
    }) }
  </NavDropdown>);
});

const GeneDropdown = withRouter(function(props) {
  const { location } = props,
    { patient_id, gene } = parse_route(location.pathname),
    alleles = patient_v_pairs[patient_id][gene];
  return (<NavDropdown title='Gene' id='gene'>
    {alleles.map(allele=> {
      const route = get_route(patient_id, gene, allele);
      return (<LinkContainer to={route} key={allele}>
        <MenuItem>
          {allele}
        </MenuItem>
      </LinkContainer>);
    }) }
  </NavDropdown>);
});

function App() {
  return(<Router>
    <div>
      <Navbar onSelect={key=>this.onSelect(key)} fixedTop>
        <Navbar.Header>
          <Navbar.Brand>
            ACME Immunology Visualization
          </Navbar.Brand>
        </Navbar.Header>
        <Nav>
          <LinkContainer exact to="/">
            <NavItem>
              Overview
            </NavItem>
          </LinkContainer>
          <PatientDropdown />
          <FamilyDropdown />
          <GeneDropdown />
        </Nav>
      </Navbar>
      <Grid fluid>
        <Route path="/:patient_id&:allele" component={BCellInterface} />
        <Route exact path="/" component={Overview} />
      </Grid>
    </div>
  </Router>);
}

ReactDOM.render(
  <App />,
  document.body.appendChild(document.createElement('div'))
)
