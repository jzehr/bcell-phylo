import React, { Component } from 'react';
import { Row, Col, Grid } from 'react-bootstrap';
import 'phylotree/phylotree.css';
const _ = require('underscore');

const d3 = require("d3");

const $ = require("jquery");

import {
  BaseAlignment,
  fastaParser,
  ScrollBroadcaster,
  SequenceAxis,
  SiteAxis
} from 'alignment.js';


require("phylotree");

const nucleotide_colors = {
    a: "LightPink",
    g: "LemonChiffon",
    t: "LightBlue",
    c: "MediumPurple",
    n: "DarkGrey",
    A: "LightPink",
    G: "LemonChiffon",
    T: "LightBlue",
    C: "MediumPurple",
    N: "DarkGrey",
    "-": "lightgrey"
  },
  nucleotide_text_colors = {
    a: "Red",
    g: "GoldenRod",
    t: "Blue",
    c: "DarkMagenta",
    n: "black",
    A: "Red",
    G: "GoldenRod",
    T: "Blue",
    C: "DarkMagenta",
    N: "black",
    "-": "DarkGrey"
  };


const colors = ['red', 'pink', 'blue', 'lightblue', 'purple', 'plum'];

function get_name(node, key) {
  key = key || 'name';
  return node[key].split('_')[0];
}

function get_size(node, key) {
  key = key || 'name';
  return +node[key].split('_')[2].split('-')[1];
}

function get_time(node, key) {
  key = key || 'name';
  return +node[key].split('_')[1].split('-')[1]-1;
}

function get_cdr3(node, key) {
  key = key || 'name';
  return node[key].split('_')[3];
}

function get_tag(node, key) {
  key = key || 'name';
  return node[key].split('_')[3];
}

function get_signature(node, key) {
  key = key || 'name';
  return node[key].split('_')[4];
}

function LegendItem(props) {
  key = key || 'name';
  const { text, color } = props;
}

function edgeStyler(element, edge) {
  if(!d3.layout.phylotree.is_leafnode(edge.target)) return null;
  const time = get_time(edge.target),
    stroke = colors[time];
  element.style("stroke", stroke);
}

function guideEdgeStyler(element, edge) {
  if(!d3.layout.phylotree.is_leafnode(edge.target)) {
    return null
  };
  const time = get_time(edge.target),
    stroke = colors[time];
  element.style("stroke", stroke);
}

class BCellPhylo extends Component {
  constructor(props) {
    super(props);

    this.column_sizes = [700, 700, 200, 200, 700];
    this.row_sizes = [20, 20, 700];
  }
  componentWillUpdate(nextProps) {
    if (nextProps.json) {
      const { site_size } = nextProps;
      this.sequence_data = fastaParser(nextProps.json.fasta)
        .map(record => {
          record.old_header = record.header;
          if (record.header.indexOf('Germline_') > -1) return record;
          record.size = get_size(record, 'header');
          record.time = get_time(record, 'header');
          record.cdr3 = get_cdr3(record, 'header');
          record.header = get_signature(record, 'header');
          return record;
        });
      const number_of_sequences = this.sequence_data.length-1;
      this.tree_size = number_of_sequences * site_size;
      this.main_tree = d3.layout
        .phylotree()
        .options({
          "left-right-spacing": "fit-to-size",
          "top-bottom-spacing": "fit-to-size",
          "show-scale": false,
          "align-tips": true,
          "show-labels": false,
          selectable: false
        })
        .style_edges(edgeStyler)
        .size([this.tree_size, this.tree_size])
        .node_circle_size(0);
      this.parsed = d3.layout.newick_parser(nextProps.json.newick);
      this.main_tree(this.parsed);

      var i = 0;
      this.main_tree.traverse_and_compute(function(n) {
        var d = 1;
        if (!n.name) {
          n.name = "Node" + i++;
        }
        if (n.children && n.children.length) {
          d += d3.max(n.children, function(d) {
            return d["count_depth"];
          });
        }
        n["count_depth"] = d;
      });

      this.main_tree.resort_children(function(a, b) {
        return a["count_depth"] - b["count_depth"];
      }, null, null, true);

      const ordered_leaf_names = this.main_tree
        .get_nodes(true)
        .filter(d3.layout.phylotree.is_leafnode)
        .map(d => d.name.split('_')[0]);

      this.sequence_data.sort((a, b) => {
        if(a.old_header.indexOf('Germline_') > -1) return -1;
        if(b.old_header.indexOf('Germline_') > -1) return 1;
      
        const a_header = a.old_header.split('_')[0],
          b_header = b.old_header.split('_')[0],
          a_index = ordered_leaf_names.indexOf(a_header),
          b_index = ordered_leaf_names.indexOf(b_header);
        return a_index - b_index;
      });
    }
  }
  componentDidUpdate() {
    if(this.props.json) {
      const { sequence_data } = this;
      const max_size = d3.max(sequence_data.slice(1).map(d=>d.size))

      this.main_tree.svg(d3.select("#alignmentjs-largeTreeAlignment")).layout();

      const guide_height = this.row_sizes[2],
        guide_width = this.column_sizes[0];

      this.guide_tree = d3.layout
        .phylotree()
        .svg(d3.select("#alignmentjs-guideTree"))
        .options({
          "left-right-spacing": "fit-to-size",
          "top-bottom-spacing": "fit-to-size",
          collapsible: false,
          transitions: false,
          "show-scale": false,
          brush: false,
          selectable: false,
          "show-labels": false,
        })
        .style_edges(guideEdgeStyler)
        .size([guide_height, guide_width])
        .node_circle_size(0);

      this.guide_tree(this.parsed).layout();

      
      const alignment_axis_width = this.column_sizes[4],
        alignment_axis_height  = this.row_sizes[0];

      d3.select("#alignmentjs-axis-div")
        .style("width", alignment_axis_width + "px")
        .style("height", alignment_axis_height + "px");

      const number_of_sites = this.sequence_data[0].seq.length,
        number_of_sequences = this.sequence_data.length,
        { site_size } = this.props,
        alignment_width = site_size * number_of_sites,
        alignment_height = site_size * number_of_sequences;

      var alignment_axis_scale = d3.scale.linear()
        .domain([1, number_of_sites])
        .range([site_size / 2, alignment_width - site_size / 2]);

      var alignment_axis_svg = d3.select("#alignmentjs-alignment-axis");
      alignment_axis_svg.html("");
      alignment_axis_svg.attr("width", alignment_width)
        .attr("height", alignment_axis_height);

      var alignment_axis = d3.svg.axis()
        .orient("top")
        .scale(alignment_axis_scale)
        .tickValues(d3.range(1, number_of_sites, 2));

      alignment_axis_svg
        .append("g")
        .attr("class", "axis")
        .attr("transform", `translate(0, ${alignment_axis_height - 1})`)
        .call(alignment_axis);

      const bar_width = this.column_sizes[3],
        bar_height = this.row_sizes[1];

      var bar_axis_svg = d3.select("#alignmentjs-bar-axis");
      bar_axis_svg.attr("width", bar_width)
        .attr("height", bar_height);

      var bar_scale = d3.scale.linear()
        .domain([0, max_size])
        .range([10, bar_width-10]);

      var bar_axis = d3.svg.axis()
        .orient("top")
        .scale(bar_scale)
        .ticks(5);

      bar_axis_svg.append("g")
        .attr("class", "axis")
        .attr("transform", `translate(0, ${alignment_axis_height - 1})`)
        .call(bar_axis);

      var bar_svg = d3.select("#alignmentjs-bar");
      bar_svg.attr("width", this.column_sizes[3])
        .attr("height", alignment_height);

      bar_svg.selectAll('rect')
        .data(sequence_data.slice(1))
        .enter()
          .append('rect')
          .attr('x', 10)
          .attr('y', function(d,i) { return i*site_size; })
          .attr('width', function(d) { return bar_scale(d.size); })
          .attr('height', site_size)
          .attr('fill', function(d) { return colors[d.time]; });
        
      const alignment_viewport_height = this.row_sizes[2],
        alignment_viewport_width = this.column_sizes[4],
        full_pixel_width = site_size * number_of_sites,
        full_pixel_height = site_size * number_of_sequences;
      const scroll_broadcaster = new ScrollBroadcaster(
        { width: full_pixel_width, height: full_pixel_height },
        { width: alignment_viewport_width, height: alignment_viewport_height },
        { x_pixel: 0, y_pixel: 0 },
        [
          "alignmentjs-axis-div",
          "germline-alignment",
          "alignmentjs-guideTree-div",
          "alignmentjs-largeTreeAlignment-div",
          "alignmentjs-labels-div",
          "alignmentjs-bar-div",
          "alignmentjs-alignment"
        ]
      );

      const guide_x_scale = d3.scale
        .linear()
        .domain([0, full_pixel_height])
        .range([0, guide_width]);
      const guide_y_scale = d3.scale
        .linear()
        .domain([0, full_pixel_height])
        .range([0, guide_height]);
      const rect = d3
        .select("#alignmentjs-guideTree")
        .append("rect")
        .attr("x", 0)
        .attr("y", 0)
        .attr("id", "guide-rect")
        .style("opacity", 0.6)
        .style("stroke-width", "1px")
        .style("stroke", "GoldenRod")
        .style("fill", "yellow")
        .attr("width", guide_x_scale(guide_width))
        .attr("height", guide_y_scale(guide_height));

      scroll_broadcaster.setListeners();

      $("#alignmentjs-guideTree-div").off("wheel");
      $("#alignmentjs-guideTree-div").on("wheel", function(e) {
        e.preventDefault();
        const guide_x = +d3.select("#guide-rect").attr("x");
        const new_guide_x = Math.min(
          Math.max(guide_x + guide_x_scale(e.originalEvent.deltaX), 0),
          guide_width - guide_x_scale(guide_width)
        );

        rect.attr("x", new_guide_x);

        const guide_y = +d3.select("#guide-rect").attr("y");
        const new_guide_y = Math.min(
          Math.max(guide_y + guide_y_scale(e.originalEvent.deltaY), 0),
          guide_height - guide_y_scale(guide_height)
        );
        rect.attr("y", new_guide_y);

        const new_x_pixel = guide_x_scale.invert(new_guide_x),
          new_y_pixel = guide_y_scale.invert(new_guide_y);
        $("#alignmentjs-largeTreeAlignment-div").scrollLeft(new_x_pixel);
        $("#alignmentjs-largeTreeAlignment-div").scrollTop(new_y_pixel);
        $("#alignmentjs-bar-div").scrollTop(new_y_pixel);

        const e_mock = {
          originalEvent: {
            deltaX: 0,
            deltaY: e.originalEvent.deltaY
          }
        };
        scroll_broadcaster.handleWheel(e_mock, "tree");
      });

      d3.select("#alignmentjs-guideTree").on("click", null);
      d3.select("#alignmentjs-guideTree").on("click", function() {
        const coords = d3.mouse(this),
          new_x_pixel = guide_x_scale.invert(coords[0]),
          new_y_pixel = guide_y_scale.invert(coords[1]),
          current_x_fraction = scroll_broadcaster.x_fraction,
          new_y_fraction = new_y_pixel / scroll_broadcaster.full_pixel_height;
        scroll_broadcaster.broadcast(current_x_fraction, new_y_fraction);
        rect.attr("x", coords[0]);
        rect.attr("y", coords[1]);
        $("#alignmentjs-largeTreeAlignment-div").scrollLeft(new_x_pixel);
        $("#alignmentjs-largeTreeAlignment-div").scrollTop(new_y_pixel);
        $("#alignmentjs-bar-div").scrollTop(new_y_pixel);
      });

      $("#alignmentjs-largeTreeAlignment-div").off("wheel");
      $("#alignmentjs-largeTreeAlignment-div").on("wheel", function(e) {
        const guide_x = +d3.select("#guide-rect").attr("x");
        const new_guide_x = Math.min(
          Math.max(guide_x + guide_x_scale(e.originalEvent.deltaX), 0),
          guide_width - guide_x_scale(guide_width)
        );
        rect.attr("x", new_guide_x);

        const guide_y = +d3.select("#guide-rect").attr("y");
        const new_guide_y = Math.min(
          Math.max(guide_y + guide_y_scale(e.originalEvent.deltaY), 0),
          guide_height - guide_y_scale(guide_height)
        );
        rect.attr("y", new_guide_y);

        const new_y_pixel = guide_y_scale.invert(new_guide_y);
        $("#alignmentjs-bar-div").scrollTop(new_y_pixel);
        const e_mock = {
          originalEvent: {
            deltaX: 0,
            deltaY: e.originalEvent.deltaY
          }
        };
        scroll_broadcaster.handleWheel(e_mock, "tree");
      });

      $("#alignmentjs-alignment").on("wheel", function(e) {
        e.preventDefault();
        const guide_y = +d3.select("#guide-rect").attr("y");
        const new_guide_y = Math.min(
          Math.max(guide_y + guide_y_scale(e.originalEvent.deltaY), 0),
          guide_height - guide_y_scale(guide_height)
        );
        rect.attr("y", new_guide_y);
        const new_y_pixel = guide_y_scale.invert(new_guide_y);
        $("#alignmentjs-bar-div").scrollTop(new_y_pixel);
        $("#alignmentjs-largeTreeAlignment-div").scrollTop(new_y_pixel);
        scroll_broadcaster.handleWheel(e, "alignment");
      });

    }
  }
  render(){
    if (!this.props.json) {
      return <div />;
    }
    const container_style = {
      display: "grid",
      gridTemplateColumns: this.column_sizes.join("px ") + "px",
      gridTemplateRows: this.row_sizes.join("px ") + "px"
    },
    legend = [
      { text: 'Visit 1, replicate 1', color: 'red' },
      { text: 'Visit 1, replicate 2', color: 'pink' },
      { text: 'Visit 2, replicate 1', color: 'blue' },
      { text: 'Visit 2, replicate 2', color: 'lightblue' },
      { text: 'Visit 3, replicate 1', color: 'purple' },
      { text: 'Visit 3, replicate 2', color: 'plum' }
    ],
    { CDR3 } = this.props.json,
    highlight_cdr3_color = (character, position, header) => {
      if(position >= CDR3[0] && position <= CDR3[1]) {
        return nucleotide_text_colors[character];
      }
      return nucleotide_colors[character];
    },
    highlight_cdr3_text_color = (character, position, header) => {
      if(position >= CDR3[0] && position <= CDR3[1]) {
        return nucleotide_colors[character];
      }
      return nucleotide_text_colors[character];
    };
    this.sequence_data.number_of_sites = this.sequence_data[0].seq.length;
    return (<Row>
      <Col xs={12}>
        <h4>Patient {this.props.json.patient}, gene {this.props.json.gene}, fragment {this.props.json.fragment}</h4>
        <div id='main_viz' style={container_style}>

          <div style={{gridArea: "1 / 1 / 3 / 4"}}>
            <svg width={1400} height={40}>
              {legend.map((d,i) => {
                return (<g key={d.text} transform={`translate(${20+175*i}, 0)`} >
                  <rect width={20} height={20} fill={d.color} />,
                  <text x={25} y={15}>{d.text}</text>
                </g>);
              }) } 
            </svg>
            <svg width={this.props.guide_size} height={this.props.guide_size} id='guide-tree' />
          </div>

        <div style={{textAlign: "center"}}><p>Family size</p></div>

        <SiteAxis
          width={this.column_sizes[4]}
          height={this.row_sizes[0]}
          sequence_data={this.sequence_data}
        />

        <div id="alignmentjs-bar-axis-div">
          <svg id="alignmentjs-bar-axis" />
        </div>

        <BaseAlignment
          width={this.column_sizes[4]}
          height={this.row_sizes[1]}
          sequence_data={[this.sequence_data[0]]}
          disableVerticalScrolling
          site_color={highlight_cdr3_color}
          text_color={highlight_cdr3_text_color}
          id='germline'
        />

        <div id="alignmentjs-guideTree-div">
          <svg id="alignmentjs-guideTree" />
        </div>

        <div
          id="alignmentjs-largeTreeAlignment-div"
          style={{ overflowX: "scroll", overflowY: "scroll" }}
        >
          <svg id="alignmentjs-largeTreeAlignment" />
        </div>

        <SequenceAxis
          width={this.column_sizes[2]}
          height={this.row_sizes[2]}
          sequence_data={this.sequence_data.slice(1)}
        />

        <div id="alignmentjs-bar-div" style={{overflowY: "scroll"}}>
          <svg id="alignmentjs-bar" />
        </div>

        <BaseAlignment
          width={this.column_sizes[4]}
          height={this.row_sizes[2]}
          sequence_data={this.sequence_data.slice(1)}
          site_color={highlight_cdr3_color}
          text_color={highlight_cdr3_text_color}
        />

        </div>
      </Col>
    </Row>);
  }
}

BCellPhylo.defaultProps = {
  'site_size': 20,
  'gridTemplateColumns': [300, 500, 0, 300],
  'gridTemplateRows': [1000]
}

module.exports = BCellPhylo;
