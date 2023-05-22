extern crate core;

use std::cmp::min;
use image::{DynamicImage, GenericImageView, Rgba, GenericImage};
use image::io::Reader as ImageReader;
use rand::Rng;

fn main() {
    println!("Hello, wave function collapse!");

    let img = ImageReader::open("img/Cave.png").unwrap().decode().unwrap();

    let mut library = Vec::new();

    let mut l_i = 0 as usize;

    for pixel in img.pixels(){
        if pixel.0 == 0 || pixel.1 == 0 {continue;}
        if pixel.0 == 15 || pixel.1 == 15 {continue;}

        let x = pixel.0;
        let y = pixel.1;

        let state = [
            img.get_pixel(x-1,y-1),img.get_pixel(x,y-1),img.get_pixel(x+1,y-1),
            img.get_pixel(x-1,y),img.get_pixel(x,y),img.get_pixel(x+1,y),
            img.get_pixel(x-1,y+1),img.get_pixel(x,y+1),img.get_pixel(x+1,y+1)
        ];

        library.push((l_i,state));
        l_i += 1;
    }

    //save library as images

    for (i,state) in library.iter().enumerate() {
        let mut d_image = DynamicImage::new_rgba8(3,3);

        d_image.put_pixel(0,0,state.1[0]);
        d_image.put_pixel(1,0,state.1[1]);
        d_image.put_pixel(2,0,state.1[2]);

        d_image.put_pixel(0,1,state.1[3]);
        d_image.put_pixel(1,1,state.1[4]);
        d_image.put_pixel(2,1,state.1[5]);

        d_image.put_pixel(0,2,state.1[6]);
        d_image.put_pixel(1,2,state.1[7]);
        d_image.put_pixel(2,2,state.1[8]);

        let p = format!("./img/lib/{}.png",i);
        d_image.save(p).unwrap();
    }

    let height = 50;
    let width = 50;
    let size = (height * width) as usize;

    let mut out = vec![Vec::with_capacity(width);height];

    for x in 0..width as u32 {
        for y in 0..height as u32 {
            out[y as usize].push(Cell::new(library.clone(),x,y));
        }
    }

    let mut entropys = vec![usize::MAX;size];
    println!("e_l:{}",entropys.len());
    let mut updates = Vec::with_capacity(size);
    entropys = entropys.iter().enumerate().map(|(i,_)|
        {
            let v_c = &out[i/width];
            let c = &v_c[i % width];
            c.collapse_coef()
        }).collect();
    let mut max_entropy = *entropys.iter().filter(|&&e| e < usize::MAX).max().unwrap();

    while max_entropy > 1{
        //println!("o_h:{} o_w:{}",out.len(),out[0].len());
        entropys = entropys.iter().enumerate().map(|(i,_)|
            {
                let y = i/width;
                let x = i%width;

                //print!("x:{} y:{} ",x,y);

                let v_c = &out[y];
                let c = &v_c[x];
                c.collapse_coef()
            }).collect();

        //println!("e_x:{}",entropys.iter().filter(|&&e| e == usize::MAX).count());

        /*
        println!("");
        entropys.iter().enumerate().for_each(|(i,&e)|{
            if i%width == 0 {println!("");}
            if e == usize::MAX {print!("X ");}
            else {
                let e_u = e as u8 + 40;
                print!("{} ",e_u as char);
            }
        });
        */
        //println!("");
        //println!("e_l:{} e:{:?}",entropys.len(),entropys);
        let (min_index,_min_value) = entropys.iter().enumerate().min_by_key(|(_,&t)|t).unwrap();

        println!("e_i:{}",min_index);

        let min_x = min_index%width;
        let min_y = min_index/width;

        out[min_y][min_x].collapse();
        updates.push((min_y,min_x));

        while updates.len() > 0 {
            let _pre_len = updates.len();

            let (u_y, u_x) = updates.pop().unwrap();
            let u_c = &out[u_y][u_x];
            let states = u_c.possible_states();

            //println!("working on: {} {} of {}",u_x,u_y,pre_len);


            let y_1 = u_y.checked_sub(1).unwrap_or(0);

            let y_2 = u_y as usize + 2;
            let y_2 = min(y_2,height);

            for v_y in &mut out[y_1..y_2] {
                let x_1 = (u_x).checked_sub(1).unwrap_or(0);

                let x_2 = u_x as usize + 2;
                let x_2 = min(x_2,width);

                for c in &mut v_y[x_1..x_2] {
                    let x_off = u_x as i32 - c.x as i32;
                    let y_off = u_y as i32 - c.y as i32;

                    if x_off == 0 && y_off == 0 {continue;}

                    //println!("x:{:02} y:{:02}",x_off,y_off);

                    let _n_s = (y_2 - y_1) * (x_2 - x_1);
                    //println!("{}",n_s);
                    //println!("{} {}",c.y,c.x);
                    if c.propograte(x_off, y_off, states.clone()) {
                        if !updates.contains(&(c.y as usize, c.x as usize)) {
                            updates.push((c.y as usize, c.x as usize));
                        }
                        //println!("wow {} {}",c.y,c.x);
                    }
                }
            }
            //println!("After: {}",updates.len());
        }
        max_entropy = *entropys.iter().filter(|&&e| e < usize::MAX).max().unwrap();
    }
    println!("done?")
}
#[derive(Clone,Debug)]
struct Cell{
    super_states: Vec<(usize,[Rgba<u8>;9])>,
    x: u32,
    y: u32,
    collapsed: bool,
}

impl Cell{
    fn new(lib : Vec<(usize,[Rgba<u8>;9])>,x : u32,y: u32) -> Cell{
        Cell{super_states: lib,x,y,collapsed:false}
    }

    fn collapse_coef(&self) -> usize{
        let s = self.super_states.len();

        if self.collapsed {
            return usize::MAX;
        }

        s
    }

    fn collapse(&mut self){
        //println!("c_l:{}",self.super_states.len());
        if self.super_states.len() == 1 {return;}

        let mut r = rand::thread_rng();
        for _ in 0..(self.super_states.len()-1) {
            let i = r.gen_range(0..self.super_states.len());
            self.super_states.remove(i);
        }

        self.collapsed = true;
        //println!("c_l_n:{}",self.super_states.len());
    }

    fn possible_states(&self) -> Vec<Rgba<u8>>{
        self.super_states.iter().map(|&(_,s)|s[4].clone()).collect()
    }

    fn propograte(&mut self,x_off:i32,y_off:i32,states: Vec<Rgba<u8>>) -> bool{
        if self.collapsed {return false;}
        let i = x_off + (y_off * 3) + 4;

        //println!("i:{} x:{} y:{}",i,x_off,y_off);

        let old_states = self.super_states.clone();
        let squares : Vec<(usize,(usize,Rgba<u8>))> = self.super_states.iter().enumerate().map(|(j,&(k,s))| (j,(k,s[i as usize]))).collect();
        let squares_filter : Vec<&usize> = squares.iter().filter(|&&(_,(_,s))| states.contains(&s)).map(|(i,_)| i).collect();


        if squares_filter.len() == old_states.len() {return false;}
        println!("x:{} y{}",self.x,self.y);
        //println!("old: l:{} {:?}",old_states.len(),old_states.iter().map(|(i,_)|i).collect::<Vec<&usize>>());
        //println!("other: l_s:{} l_f:{}",squares.len(),squares_filter.len());

        //println!("s:{} s_f:{} s_t:{}",squares.len(),squares_filter.len(),states.len());
        println!("i: {} : {:?}",i,states);
        //println!("squares:{:?}\nstates:{:?}\noff:{},{}\nfilter:{:?}",self.super_states,states,x_off,y_off,squares_filter);

        let indexes : Vec<(usize,&(usize,[Rgba<u8>;9]))> = old_states.iter().enumerate().collect();
        let removed : Vec<(&usize,&usize)> = indexes.iter().filter(|(i,_)|!squares_filter.contains(&i)).map(|(i,(r,_))|(i,r)).collect();

        for j in 0..196 {
            if j % 14 == 0
            {
                println!("");
                print!("# ");
            }
            let r_i : Vec<&(&usize,&usize)> = removed.iter().filter(|(_,&k)|j == k).collect();
            let s_i : Vec<&&usize> = squares_filter.iter().filter(|&&&k| j == k).collect();
            if r_i.len() > 0 {print!("x")}
            else if s_i.len() > 0{print!("o")}
            else { print!(".") }
        }
        println!("");
        //println!("keep:{:?}",squares_filter);
        if squares_filter.len() == 0 {
            panic!("Will remove all states!\n
            x: {} y:{}\n
            squares:{:?}\n
            states:{:?}\n
            filter:{:?}\n
            off set x:{} y:{}\n
            removed:{:?}
            ",self.x,self.y,squares,states,squares_filter,x_off,y_off,removed);}
        //println!("o:{}",self.super_states.len());

        //self.super_states = self.super_states.iter().enumerate().filter(|(i,_)| squares_filter.contains(&i)).map(|(_,&s_s)| s_s).collect();

        for (r,_) in removed.iter().rev() {
            self.super_states.swap_remove(**r);
        }

        if self.super_states.len() == 0 {panic!("No States??");}

        //println!("n:{}",self.super_states.len());
        //println!("c:{}",changed);

        true

    }
}


