use std::env;

fn print_example() {
    println!("Example usage: ");
    println!("cargo run imgx imgy scale");
    println!("cargo run 100 100 1.0");
    panic!()
}

struct Settings {
    imgx: u32,
    imgy: u32,
    scale: f32,
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
        let outname = args[4].clone();
      
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

fn main() {

    let settings = Settings::new(env::args().collect());

    // Create a new ImgBuf with width: imgx and height: imgy
    let mut imgbuf = image::ImageBuffer::new(settings.imgx, settings.imgy);

    // Iterate over the coordinates and pixels of the image
    for (i, j, pixel) in imgbuf.enumerate_pixels_mut() {
        let x = settings.scale * ((i as f32 + 0.5) / settings.imgx as f32 - 0.5);
        let y = settings.scale * ((j as f32 + 0.5) / settings.imgy as f32 - 0.5);

        //*pixel = image::Rgb([x as u8 ,y as u8, 0]);
        if x * x + y * y < 0.25 {
            *pixel = image::Rgb(to_color_rgb(x.abs(),y.abs(),0.0));
        }
    }

    // Save the image
    imgbuf.save(settings.outname).unwrap();
}
