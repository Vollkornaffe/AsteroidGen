use cgmath::vec2;
use cgmath::vec3;
use cgmath::InnerSpace;
use rand::prelude::*;
use std::env;

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
    [(255.0 * r) as u8, (255.0 * g) as u8, (255.0 * b) as u8]
}

struct Node {
    pos: cgmath::Vector2<f32>,
    size: f32,
}

#[derive(Clone, Copy)]
struct Edge {
    a: usize,
    b: usize,
    length: f32,
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

    fn add_node(&mut self, pos: cgmath::Vector2<f32>, size: f32) {
        self.nodes.push(Node { pos, size });
    }

    fn add_possible_edge(&mut self, a: usize, b: usize, size: f32) {
        if self.nodes.len() <= a || self.nodes.len() <= b {
            panic!(
                "Invalid a or b for graph! #nodes: {}, a: {}, b: {}",
                self.nodes.len(),
                a,
                b
            );
        }

        self.possible_edges.push(Edge {
            a,
            b,
            length: (self.nodes[a].pos - self.nodes[b].pos).magnitude(),
            size,
        });
    }

    fn dist_to_nodes(&self, q: cgmath::Vector2<f32>) -> f32 {
        let mut min_dist = std::f32::MAX;
        for node in &self.nodes {
            let dist = (q - node.pos).magnitude() - node.size;
            if dist < min_dist {
                min_dist = dist;
            }
        }
        min_dist
    }

    fn dist_to_edges(&self, q: cgmath::Vector2<f32>) -> f32 {
        let mut min_dist = std::f32::MAX;
        for edge in &self.edges {
            let a = self.nodes[edge.a].pos;
            let b = self.nodes[edge.b].pos;
            let qa = q - a;
            let ba = b - a;
            let h = num::clamp((qa.dot(ba)) / (ba.dot(ba)), 0.0, 1.0);
            let dist = (qa - h * ba).magnitude();
            if dist < min_dist {
                min_dist = dist;
            }
        }
        min_dist
    }
}

struct Triangle {
    a: cgmath::Vector2<f32>,
    b: cgmath::Vector2<f32>,
    c: cgmath::Vector2<f32>,
    det: f32,
}

impl Triangle {
    fn new(a: cgmath::Vector2<f32>, b: cgmath::Vector2<f32>, c: cgmath::Vector2<f32>) -> Self {
        Self {
            a,
            b,
            c,
            det: (b.y - c.y) * (a.x - c.x) + (c.x - b.x) * (a.y - c.y),
        }
    }

    fn contains(&self, p: cgmath::Vector2<f32>) -> bool {
        let bary_coords = self.cart_to_bary(p);
        bary_coords.x < 0.0 || bary_coords.y < 0.0 || bary_coords.z < 0.0
    }

    fn bary_to_cart(&self, coords: cgmath::Vector3<f32>) -> cgmath::Vector2<f32> {
        self.a * coords.x + self.b * coords.y + self.c * coords.z
    }

    fn cart_to_bary(&self, coords: cgmath::Vector2<f32>) -> cgmath::Vector3<f32> {
        let mu_1 = ((self.b.y - self.c.y) * (coords.x - self.c.x)
            + (self.c.x - self.b.x) * (coords.y - self.c.y))
            / self.det;
        let mu_2 = ((self.c.y - self.a.y) * (coords.x - self.c.x)
            + (self.a.x - self.c.x) * (coords.y - self.c.y))
            / self.det;
        let mu_3 = 1.0 - mu_1 - mu_2;
        vec3(mu_1, mu_2, mu_3)
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
        graph.add_node(vec2(x, y), 0.01);

        if graph.nodes.len() > 7 || (exit < 0.1 && graph.nodes.len() > 3) {
            break;
        }
    }

    // calculate dulaunay triangulation so we can have a planar graph
    let delaunay_points: Vec<delaunator::Point> = graph
        .nodes
        .iter()
        .map(|n| delaunator::Point {
            x: n.pos.x as f64,
            y: n.pos.y as f64,
        })
        .collect();
    let delaunay_triangulation =
        delaunator::triangulate(&delaunay_points).expect("Delaunator failed.");

    println!("{:?}", delaunay_triangulation.hull);

    let mut triangles = Vec::new();

    // then use the triangulation to generate some possible edges
    for i in 0..(delaunay_triangulation.triangles.len() / 3) {
        let a_i = delaunay_triangulation.triangles[i * 3 + 0];
        let b_i = delaunay_triangulation.triangles[i * 3 + 1];
        let c_i = delaunay_triangulation.triangles[i * 3 + 2];

        let a = graph.nodes[a_i].pos;
        let b = graph.nodes[b_i].pos;
        let c = graph.nodes[c_i].pos;

        triangles.push(Triangle::new(a, b, c));

        graph.add_possible_edge(
            delaunay_triangulation.triangles[i * 3 + 0],
            delaunay_triangulation.triangles[i * 3 + 1],
            0.005,
        );
        graph.add_possible_edge(
            delaunay_triangulation.triangles[i * 3 + 1],
            delaunay_triangulation.triangles[i * 3 + 2],
            0.005,
        );
        graph.add_possible_edge(
            delaunay_triangulation.triangles[i * 3 + 2],
            delaunay_triangulation.triangles[i * 3 + 0],
            0.005,
        );
    }

    // sort the edges by length so that the tree is minimal
    graph
        .possible_edges
        .sort_by(|a, b| b.length.partial_cmp(&a.length).unwrap());

    // create MST from triangulation
    // need to init with one edge
    graph.edges.push(graph.possible_edges[0]);
    // track which nodes are in the MST so far
    let mut connected = vec![graph.edges[0].a, graph.edges[0].b];
    // just add all edges that don't create a loop.
    loop {
        let mut smt_new = false;
        for &possible_edge in &graph.possible_edges {
            let a = possible_edge.a;
            let b = possible_edge.b;
            match (connected.contains(&a), connected.contains(&b)) {
                (true, true) => continue,
                (false, true) => {}
                (true, false) => {}
                (false, false) => continue,
            }
            smt_new = true;
            connected.push(possible_edge.a);
            connected.push(possible_edge.b);
            graph.edges.push(possible_edge);
            continue;
        }
        if !smt_new {
            break;
        }
    }

    // Create a new ImgBuf with width: imgx and height: imgy
    let mut imgbuf = image::ImageBuffer::new(settings.imgx, settings.imgy);

    // Iterate over the coordinates and pixels of the image
    for (i, j, pixel) in imgbuf.enumerate_pixels_mut() {
        let x = (i as f32 + 0.5) / settings.imgx as f32 - 0.5;
        let y = (j as f32 + 0.5) / settings.imgy as f32 - 0.5;

        let pixel_pos = vec2(x, y);

        if graph.dist_to_edges(pixel_pos) < 0.0 {
            *pixel = image::Rgb(to_color_rgb(0.0, 0.0, 1.0));
        }
        if graph.dist_to_nodes(pixel_pos) < 0.0 {
            *pixel = image::Rgb(to_color_rgb(1.0, 0.0, 0.0));
        }
    }

    // Save the image
    imgbuf.save(settings.outname).unwrap();
}
