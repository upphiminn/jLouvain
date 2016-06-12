/* 
 Example of jLouvain used from node.js

 */
var jLouvain = require('../jLouvain')

var node_data = ['id1', 'id2', 'id3', 'id4']
	, edge_data = [
		{source: 'id1', target:'id2', weight: 10.0}
		, {source: 'id2', target:'id3', weight: 20.0}
		, {source: 'id3', target:'id1', weight: 30.0}
	]
	, init_part = {'id1':0, 'id2':0, 'id3': 1}

var community = jLouvain().nodes(node_data).edges(edge_data).partition_init(init_part)
	, result  = community()
	
console.log('result is:', result)
