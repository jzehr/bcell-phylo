import React, { Component } from 'react';
import ReactTable from "react-table";
import { withRouter } from "react-router-dom";
const d3 = require("d3");
import 'react-table/react-table.css'

const PATIENT_IDS = ["28729", "48689", "67029", "77612", "78202", "93954", "99361", "99682", "GJS"];

const ExpansionTable = withRouter(function(props){
  const { history } = props;
  return (<ReactTable
    data={props.data}
    defaultSortDesc={true}
    className='-highlight'
    columns={[
      { Header: "Patient ID", accessor: 'patient_id', headerStyle: {fontWeight: 'bold'}},
      { Header: "V-gene allele", accessor: 'vgene_allele', headerStyle: {fontWeight: 'bold'}},
      { Header: "Visit 1 size", accessor: 'visit1', headerStyle: {fontWeight: 'bold'}},
      { Header: "Visit 2 size", accessor: 'visit2', headerStyle: {fontWeight: 'bold'}},
      { Header: "Visit 3 size", accessor: 'visit3', headerStyle: {fontWeight: 'bold'}}
    ]}
    getTrProps={(state, rowInfo, column) => {
      if (!rowInfo) return { className: '' };
      const { patient_id, vgene_allele } = rowInfo.row,
        route = `/${patient_id}&${vgene_allele}`;
      return {
        onClick: ()=>history.push(route)
      }
    }}
  />);
});

class Overview extends Component {
  constructor(props) {
    super(props);
    this.state = {
      data: null,
      patient_id: 'all',
      vgene: 'all'
    };
  }
  componentDidMount() {
    d3.csv('data/cluster_information.csv', (error, data) => {
      data.forEach(row=>{
        row.visit1 = +row.visit1;
        row.visit2 = +row.visit2;
        row.visit3 = +row.visit3;
      });
      this.setState({data: data});
    });
  }
  render() {
    if (!this.state.data) return <div/>;
    const data = this.state.data
      .filter(datum => {
        const desired_patient_id = this.state.patient_id == 'all' ?
          true :
          datum.patient_id == this.state.patient_id,
        desired_vgene = this.state.vgene == 'all' ?
          true :
          +datum.vgene_allele[1] == this.state.vgene,
        is_desired = desired_patient_id && desired_vgene;
        return is_desired;
      });
    return (<div style={{display: "flex", justifyContent: "center"}}>
      <div>
        <div style={{display: "flex", justifyContent: "space-around"}}> 
          <b> FILTERS - </b>
          <label>
            Patient ID: 
            <select value={this.state.patient_id} onChange={e=>this.setState({patient_id: e.target.value})}>
              <option value='all'>All</option>
              {PATIENT_IDS.map(patient_id => <option value={patient_id}>{patient_id}</option>)}
            </select>
          </label>
          <label>
            V-gene: 
            <select value={this.state.vgene} onChange={e=>this.setState({vgene: e.target.value})}>
              <option value='all'>All</option>
              {[1,2,3,4,5,6].map(vgene => <option value={vgene}>{vgene}</option>)}
            </select>
          </label>
        </div>
        <ExpansionTable data={data} /> 
      </div>
      <div style={{padding: 20}}>
        <h2>BCell-Phylo</h2>
        <p>Phylogenetic visualization of immune repertoires.</p>
        <h3>Overview</h3>
        <p>The present page.</p>
        <p><b>Click</b> a row of the table to view results for that patient/V-gene combination.</p>
        <p><b>Filter</b> a column by selecting the appropriate filter from the dropdown at the top.</p>
        <p><b>Click</b> a column header to sort a given column in descending order.</p>
        <p><b>Size</b> denotes largest clades found which are composed of 90% or more from a particular visit (to detect clonal expansion).</p>
        <h3>Visualization</h3>
        <p>Alignments, phylogenetic trees, and annotated regions are displayed.</p>
      </div>
    </div>);
  }
}

module.exports = Overview;

