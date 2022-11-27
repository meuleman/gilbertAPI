import React from "react";
import * as d3 from "d3"


class BarD3 extends React.Component {

    render() {
        var data = this.props.barProps.regionalVals;
        // console.log(data)
        
        if (data.length !== 0) {
            var margin = {top: 0, right: 50, left: 50, bottom: 10},
                height = 1000,
                width = 800;
            
            
            
            d3.select("svg").remove()
            
            if (document.getElementById("barChart").getElementsByTagName('*').length !== 0) {
                document.getElementById("barChart").removeChild(document.getElementById("barChart").lastElementChild);
            }
            //
            // d3.select("div").remove()
            var barContainer = d3.select("#barChart");
            

            var barsvg = barContainer
                        .append("svg")
                            .attr("width", width)
                            .attr("height", height)
                        .append("g")
                            .attr("transform", "translate(" + margin.left + "," + margin.top + ")");

            var metrics = Object.keys(data[0]).slice(0, 3)
            //var regions = d3.map(data, function(d){return(d.region)}).keys()
            var regionSums = data.map(x => {
                return Object.values(x).slice(0, 3).reduce(
                    (regionTotal, metricVal) => regionTotal + metricVal, 0
                );
            });

            var regionScores = data.map(x => {
                return 1 - Object.values(x).slice(0, 3).reduce(
                    (scoreIterator, metricVal) => scoreIterator * (1 - metricVal), 1
                );
            });

            var scaledData = data.map((x, r) => {
                var regionScaled = {}
                metrics.map((m) => {
                    regionScaled[m] = x[m] / regionSums[r] * regionScores[r]
                });
                regionScaled["region"] = r;
                return regionScaled
            })
            
            // X axis
            var x = d3.scaleLinear()
                .domain([0, Math.max(...regionSums)])
                .range([0, width - margin.left - margin.right]);
            barsvg.append("g")
                .style("font", "8px times")
                .attr("transform", `translate(0, ${height})`)
                .call(d3.axisBottom(x));

            // Y axis
            var y = d3.scaleBand()
                .domain(data.map(d => d.region))
                .range([0, height - margin.top - margin.bottom])
                .padding(0.1);
            barsvg.append("g")
                .style("font", "8px times")
                // .attr("transform", "translate(0," + width + ")")
                .call(d3.axisLeft(y));

            var color = d3.scaleOrdinal()
                .domain(metrics)
                .range(['#a6cee3','#1f78b4','#b2df8a'])

            var stackedData = d3.stack().keys(metrics)(scaledData)

            
            var tooltip = barContainer
                .append("div")
                .style("opacity", 0)
                .attr("class", "tooltip")
                .style("background-color", "white")
                .style("border", "solid")
                .style("border-width", "1px")
                .style("border-radius", "5px")
                .style("padding", "10px")
                .style("position", "absolute")
                .style("display", "inline")
                .style("z-index", 1)
                

            var mouseover = function(d, i) {
                var metricName = d3.select(this.parentNode).datum().key;
                var regionNum = i.data.region;
                tooltip
                    .html("Region: " + regionNum + "<br>Metric: " + metricName)
                    .style("opacity", 1);
            }

            var mousemove = function(e) {
                tooltip.style("left", (e.pageX+20) + "px");
                tooltip.style("top", (e.pageY) + "px");
                // console.log(tooltip.style.top, tooltip.style.left)
            }
            var mouseleave = function(d) {
                tooltip
                    .style("opacity", 0)
            }

            barsvg.append("g")
                .selectAll("g")
                .data(stackedData)
                .enter().append("g")
                    .attr("fill", function(d) {return color(d.key); })
                    .selectAll("rect")
                    .data(function(d) {
                        return d;
                    })
                    .enter().append("rect")
                        .attr("x", function(d) {
                            return x(d[0]); 
                        })
                        .attr("y", function(d) {
                            // console.log(y(d.data.region))
                            return y(d.data.region); 
                        })
                        .attr("width", function(d) {
                            // console.log(x(d[1]) - x(d[0]))
                            return x(d[1]) - x(d[0]);
                        })
                        .attr("height", y.bandwidth)
                    .on("mouseover", mouseover)
                    .on("mousemove", mousemove)
                    .on("mouseleave", mouseleave)


        }
        
        return (
            
            <div id="barChart">

                {/*<div>
                    {document.getElementById('barChart')}
                </div>*/}
            </div>
            
            
        )
    }
};

export default BarD3
