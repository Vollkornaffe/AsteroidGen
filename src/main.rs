use std::env;
use rand::prelude::*;

fn print_example() {
    println!("Example usage: ");
    println!("cargo run imgx imgy scale node_size edge_size");
    println!("cargo run 100 100 0.2 0.1 1.0");
    panic!()
}

struct Settings {
    imgx: u32,
    imgy: u32,
    outname: String,
}

impl Settings {
    fn new(args: Vec<String>) -> Self {

        if args.len() < 4 {
            println!("Not enough arguments given.");
            print_example();
        }

        let imgx = args[1].parse::<u32>().expect("imgx couldn't be parsed");
        let imgy = args[2].parse::<u32>().expect("imgy couldn't be parsed");
        let outname = args[3].clone();
      
        Self {
            imgx,
            imgy,
            outname,
        }
    }
}

fn to_color_rgb(r: f32, g: f32, b: f32) -> [u8; 3] {
    [
        (255.0 * r) as u8,
        (255.0 * g) as u8,
        (255.0 * b) as u8,
    ]
}

struct Node {
    pos_x: f32,
    pos_y: f32,
    size: f32,
}

#[derive(Clone, Copy)]
struct Edge {
    a: usize,
    b: usize,
    size: f32,
}

struct Graph {
    nodes: Vec<Node>,
    possible_edges: Vec<Edge>,
    edges: Vec<Edge>,
}

impl Graph {
    fn new() -> Self {
        Self {
            nodes: Vec::new(),
			possible_edges: Vec::new(),
            edges: Vec::new(),
        }
    }

    fn add_node(&mut self, pos_x: f32, pos_y: f32, size: f32) {
        self.nodes.push(Node { pos_x, pos_y, size });
    }

    fn add_possible_edge(&mut self, a: usize, b: usize, size: f32) {
        if self.nodes.len() <= a
        || self.nodes.len() <= b {
            panic!("Invalid a or b for graph! #nodes: {}, a: {}, b: {}", self.nodes.len(), a, b);
        }

        self.possible_edges.push(Edge { a, b, size });
    }

	fn choose_edge(&mut self, i: usize) {
		self.edges.push(self.possible_edges[i]);
		self.possible_edges.remove(i);
	}

    fn dist_to_nodes(&self, q_x: f32, q_y: f32) -> f32 {
        let mut min_dist = std::f32::MAX;
        for node in &self.nodes {
            let dist_x = q_x - node.pos_x;
            let dist_y = q_y - node.pos_y;
            let dist = (dist_x*dist_x + dist_y*dist_y).sqrt() - node.size;
            if dist < min_dist {
                min_dist = dist;
            }
        }
        min_dist
    }

    fn dist_to_edges(&self, q_x: f32, q_y: f32) -> f32 {
        let mut min_dist = std::f32::MAX;
        for edge in &self.edges {
            let a_x = self.nodes[edge.a].pos_x;
            let a_y = self.nodes[edge.a].pos_y;
            let b_x = self.nodes[edge.b].pos_x;
            let b_y = self.nodes[edge.b].pos_y;
            let qa_x = q_x - a_x;
            let qa_y = q_y - a_y;
            let ba_x = b_x - a_x;
            let ba_y = b_y - a_y;
            let h = num::clamp(
                (qa_x * ba_x + qa_y * ba_y) / (ba_x * ba_x + ba_y * ba_y),
                0.0, 1.0,
            );
            let dist_x = qa_x - ba_x * h;
            let dist_y = qa_y - ba_y * h;
            let dist = (dist_x * dist_x + dist_y * dist_y).sqrt() - edge.size;
            if dist < min_dist {
                min_dist = dist;
            }
        }
        min_dist
    }
}

fn main() {

    let settings = Settings::new(env::args().collect());

    let mut graph = Graph::new();

	// add some random nodes
    loop {
        let mut rng = rand::thread_rng();
        let exit: f32 = rng.gen();

        let x = rng.gen::<f32>() - 0.5;
        let y = rng.gen::<f32>() - 0.5;
        println!("x: {}, y: {}", x, y);
        graph.add_node(x, y, 0.01);

        if exit < 0.1 && graph.nodes.len() > 3 {
            break;
        }
    }

	// calculate dulaunay triangulation so we can have a planar graph
	let delaunay_points: Vec<delaunator::Point> = graph.nodes.iter().map(
		|n| delaunator::Point {
			x: n.pos_x as f64,
			y: n.pos_y as f64,
		}
	).collect();
	let delaunay_triangulation =
		delaunator::triangulate(&delaunay_points)
		.expect("Delaunator failed.");

	println!("{:?}", delaunay_triangulation.triangles);
	println!("{:?}", delaunay_triangulation.halfedges);

	// then use the triangulation to generate some possible
	for i in 0..delaunay_triangulation.halfedges.len() {
		graph.add_possible_edge(
			delaunay_triangulation.halfedges[i],
			delaunay_triangulation.triangles[i],
			0.005,
		);
	}
	graph.edges = graph.possible_edges.clone();
	
    // Create a new ImgBuf with width: imgx and height: imgy
    let mut imgbuf = image::ImageBuffer::new(settings.imgx, settings.imgy);

    // Iterate over the coordinates and pixels of the image
    for (i, j, pixel) in imgbuf.enumerate_pixels_mut() {
        let x = (i as f32 + 0.5) / settings.imgx as f32 - 0.5;
        let y = (j as f32 + 0.5) / settings.imgy as f32 - 0.5;

        if graph.dist_to_edges(x,y) < 0.0 {
            *pixel = image::Rgb(to_color_rgb(0.0,0.0,1.0));
        }
        if graph.dist_to_nodes(x,y) < 0.0 {
            *pixel = image::Rgb(to_color_rgb(1.0,0.0,0.0));
        }
    }

    // Save the image
    imgbuf.save(settings.outname).unwrap();
}
