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
    scale: f32,
    node_size: f32,
    edge_size: f32,
    outname: String,
}

impl Settings {
    fn new(args: Vec<String>) -> Self {

        if args.len() < 3 {
            println!("Not enough arguments given.");
            print_example();
        }

        let imgx = args[1].parse::<u32>().expect("imgx couldn't be parsed");
        let imgy = args[2].parse::<u32>().expect("imgy couldn't be parsed");
        let scale = args[3].parse::<f32>().expect("scale couldn't be parsed");
        let node_size = args[4].parse::<f32>().expect("node_size couldn't be parsed");
        let edge_size = args[5].parse::<f32>().expect("edge_size couldn't be parsed");
        let outname = args[6].clone();
      
        Self {
            imgx,
            imgy,
            scale,
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
struct Edge {
    a: usize,
    b: usize,
    size: f32,
}

struct Graph {
    nodes: Vec<Node>,
    edges: Vec<Edge>,
}

impl Graph {
    fn new() -> Self {
        Self {
            nodes: Vec::new(),
            edges: Vec::new(),
        }
    }

    fn add_node(&mut self, pos_x: f32, pos_y: f32, size: f32) {
        self.nodes.push(Node { pos_x, pos_y, size });
    }

    fn add_edge(&mut self, a: usize, b: usize, size: f32) {
        if self.nodes.len() <= a
        || self.nodes.len() <= b {
            panic!("Invalid a or b for graph! #nodes: {}, a: {}, b: {}", self.nodes.len(), a, b);
        }

        self.edges(push(Edge { a, b, size }));
    }

    fn dist_to_nodes(&self, q_x: f32, q_y: f32) -> f32 {
        for node in self.nodes {
            let dist_x = q_x - node.pos_x;
            let dist_y = q_y - node.pos_y;
            let dist = (dist_x*dist_x + dist_y*dist_y).sqrt() - nodes.size;
        }
    }

    fn dist_to_edges(&self, q_x: f32, q_y: f32) -> f32 {
        for edge in self.edges {
            let a_x = self.nodes[edge.a].pos_x;
            let a_y = self.nodes[edge.a].pos_y;
            let b_x = self.nodes[edge.b].pos_x;
            let b_y = self.nodes[edge.b].pos_y;
            let d_a_x = q_x - a_x;
            let d_a_y = q_y - a_y;
            let d_b_x = q_x - b_x;
            let d_b_y = q_y - b_y;
            let h = num::clamp(
                () / (),
                0.0,
                1.0,
            );
        }

        float sdCapsule( vec3 p, vec3 a, vec3 b, float r )
        {
            vec3 pa = p - a;
            vec3 ba = b - a;
            float h = clamp( dot(pa,ba)/dot(ba,ba), 0.0, 1.0 );
            return length( pa - ba*h ) - r;
        }
        let dist_x = node.pos_x - q_x;
        let dist_y = node.pos_y - q_y;
        let dist = (dist_x*dist_x + dist_y*dist_y).sqrt() - nodes.size;

    }
}

fn main() {

    let settings = Settings::new(env::args().collect());

    let mut graph = Graph::new();

    graph.add_node(0.0, 0.0);
    graph.add_node(0.0, 0.5);
    graph.add_edge(0, 1);

    // Create a new ImgBuf with width: imgx and height: imgy
    let mut imgbuf = image::ImageBuffer::new(settings.imgx, settings.imgy);

    // Iterate over the coordinates and pixels of the image
    for (i, j, pixel) in imgbuf.enumerate_pixels_mut() {
        let x = settings.scale * ((i as f32 + 0.5) / settings.imgx as f32 - 0.5);
        let y = settings.scale * ((j as f32 + 0.5) / settings.imgy as f32 - 0.5);



        if x * x + y * y < 0.25 {
            *pixel = image::Rgb(to_color_rgb(x.abs(),y.abs(),0.0));
        }
    }

    // Save the image
    imgbuf.save(settings.outname).unwrap();
}
