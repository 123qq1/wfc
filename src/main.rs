use image::GenericImageView;
use image::io::Reader as ImageReader;

fn main() {
    println!("Hello, wave function collapse!");

    let img = ImageReader::open("img/Cave.png").unwrap().decode().unwrap();

    let mut library = Vec::new();

    for pixel in img.pixels(){
        if pixel.0 == 0 || pixel.1 == 0 {continue;}
        if pixel.0 == 15 || pixel.1 == 15 {continue;}

        let x = pixel.0;
        let y = pixel.1;

        let state = [
            img.get_pixel(x-1,y-1),img.get_pixel(x,y-1),img.get_pixel(x+1,y-1),
            img.get_pixel(x-1,y),pixel,img.get_pixel(x+1,y),
            img.get_pixel(x-1,y+1),img.get_pixel(x,y+1),img.get_pixel(x+1,y+1)
        ];

        library.push(state);

        print!("pixel {} {}, ",pixel.0,pixel.1);
    }
}
