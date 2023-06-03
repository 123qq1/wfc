use std::cmp::min;
use std::hash::Hash;
use std::iter::zip;
use image::{DynamicImage, GenericImageView, Rgba, GenericImage};
use image::io::Reader as ImageReader;
use rand::Rng;

fn main() {
    println!("Hello, wave function collapse!");

    let img = ImageReader::open("img/Knot.png").unwrap().decode().unwrap();

    let mut library = Vec::new();

    let mut l_i = 0 as usize;

    let state_size = ((img.height()-2) * (img.width()-2)) as usize;
    let req_bytes = ((state_size + 64 - 1) / 64);

    let mut n_base =
        [
            vec![0;req_bytes], vec![0;req_bytes],vec![0;req_bytes],
            vec![0;req_bytes], vec![0;req_bytes],vec![0;req_bytes],
            vec![0;req_bytes], vec![0;req_bytes],vec![0;req_bytes]
        ];

    for pixel in img.pixels(){
        if pixel.0 == 0 || pixel.1 == 0 {continue;}
        if pixel.0 == (img.width()-1) || pixel.1 == (img.height()-1) {continue;}

        let x = pixel.0;
        let y = pixel.1;

        let state = [
            img.get_pixel(x-1,y-1),img.get_pixel(x,y-1),img.get_pixel(x+1,y-1),
            img.get_pixel(x-1,y),img.get_pixel(x,y),img.get_pixel(x+1,y),
            img.get_pixel(x-1,y+1),img.get_pixel(x,y+1),img.get_pixel(x+1,y+1)
        ];

        library.push(State{state,neighbours:n_base.clone(),index:l_i,size:state_size});
        l_i += 1;
    }

    println!("Sampling done {}",library.len());

    let lib_c = library.clone();

    //fill out the library with neighbours
    for state in &mut library {
        for neighbour in &lib_c {
            for i in 0..9 {
                if i == 4 {continue;}
                state.neighbour_check(neighbour.clone(),i);
            }
        }
    }

    println!("neighbour done");

    /*

    println!("");
    for state in &library {
        print!("state {:03} :",state.index);
        for neighbour in &state.neighbours {
            print!("{:03} ",neighbour.len());
        }
        println!("");
    }

     */
    let height = 40;
    let width = 40;
    let size = (height * width) as usize;

    let mut out = vec![Vec::with_capacity(width);height];

    for x in 0..width as u32 {
        for y in 0..height as u32 {
            out[y as usize].push(Cell::new(library.clone(),x,y,library.len()));
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

        println!("e_x:{}",entropys.iter().filter(|&&e| e == usize::MAX).count());

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

        let min_x = min_index%width;
        let min_y = min_index/width;

        out[min_y][min_x].collapse();
        updates.push((min_y,min_x));

        while updates.len() > 0 {
            let _pre_len = updates.len();

            let (u_y, u_x) = updates.pop().unwrap();

            //println!("working on: {} {} of {}",u_x,u_y,pre_len);

            let u_c = out[u_y][u_x].clone();

            let y_1 = u_y.checked_sub(1).unwrap_or(0);

            let y_2 = u_y as usize + 2;
            let y_2 = min(y_2,height);

            for v_y in &mut out[y_1..y_2] {
                let x_1 = (u_x).checked_sub(1).unwrap_or(0);

                let x_2 = u_x as usize + 2;
                let x_2 = min(x_2,width);

                for c in &mut v_y[x_1..x_2] {
                    let x_off = c.x as i32 - u_x as i32;
                    let y_off = c.y as i32 - u_y as i32 ;

                    let i = x_off+(y_off*3)+4;

                    if x_off == 0 && y_off == 0 {continue;}

                    let x_c_off = c.x as i32 - min_x as i32;
                    let y_c_off = c.y as i32 - min_y as i32;

                    //if x_c_off.abs() > 1 {continue;}
                    //if y_c_off.abs() > 1 {continue;}

                    //println!("x:{:02} y:{:02}",x_off,y_off);

                    let _n_s = (y_2 - y_1) * (x_2 - x_1);
                    //println!("{}",n_s);
                    //println!("{} {}",c.y,c.x);
                    if c.propograte( u_c.possible_states(i as usize).clone()) {
                        if !updates.contains(&(c.y as usize, c.x as usize)) {
                            //println!("pushed {} {}",c.x,c.y);
                            updates.push((c.y as usize, c.x as usize));
                        }
                        //println!("wow {} {}",c.y,c.x);
                    }
                }
            }
            //println!("After: {}",updates.len());
        }
        max_entropy = *entropys.iter().filter(|&&e| e < usize::MAX).max().unwrap_or(&0);
    }
    println!("done?");

    let mut d_image = DynamicImage::new_rgba8(width as u32,height as u32);

    for w in 0..width {
        for h in 0..height {

            d_image.put_pixel(w as u32,h as u32,out[h][w].state());

        }
    }
    let p = "./img/out.png";
    d_image.save(p).unwrap();
}

#[derive(Clone,Debug)]
struct State{
    index : usize,
    neighbours: [Vec<usize>;9],
    state: [Rgba<u8>;9],
    size: usize
}

impl PartialEq for State{
    fn eq(&self, other: &Self) -> bool {

        self.state == other.state

    }
}

impl State{

    fn neighbour_check(&mut self, other :State, i : usize){
        let x_self_off = 2 - (i%3);
        let y_self_off = 2 - (i/3);

        let x_other_off = i%3;
        let y_other_off = i/3;

        let mut self_slice = Vec::new();
        let mut other_slice = Vec::new();

        for (p_i,&pix) in self.state.iter().enumerate() {
            let p_x = p_i%3;
            let p_y = p_i/3;

            if p_x != x_self_off && p_y != y_self_off {
                self_slice.push(pix.clone());
            }
        }

        for (p_i,&pix) in other.state.iter().enumerate() {
            let p_x = p_i%3;
            let p_y = p_i/3;

            if p_x != x_other_off && p_y != y_other_off {
                other_slice.push(pix.clone());
            }
        }

        //println!("self: {:?}",self_slice.len());
        //println!("other: {:?}",other_slice.len());
        //println!("off : {} {} , {} {}",x_self_off,y_self_off,x_other_off,y_other_off);

        if other_slice == self_slice
        {
            let b_i = other.index/64;
            let b_o = other.index%64;

            //while self.neighbours[i].len() <= b_i {self.neighbours[i].push(0)}

            //println!("bitting {} {} {}",i,b_i,b_o);

            self.neighbours[i][b_i] |= 1<<b_o;
        }
    }
}

#[derive(Clone,Debug)]
struct Cell{
    lib: Vec<State>,
    super_states: Vec<usize>,
    x: u32,
    y: u32,
    collapsed: bool,
    size: usize,
}

impl Cell{
    fn new(lib : Vec<State>,x : u32,y: u32,size: usize) -> Cell{

        let len = (size + 64 - 1)/64;

        let mut c = Cell{lib,x,y,collapsed:false,size,super_states:vec![usize::MAX;len]};

        c.super_states[len-1] = (1 << (size%64))-1;

        c
    }

    fn collapse_coef(&self) -> usize{
        if self.collapsed {
            return usize::MAX;
        }
        let s = self.super_states.iter().map(|u|u.count_ones()).reduce(|a,b|a+b);
        s.unwrap() as usize
    }

    fn collapse(&mut self){
        //println!("c_l:{}",self.super_states.len());

        let mut indexes = Vec::new();

        for index in 0..self.super_states.len() {
            for shift in 0..64 {
                let bit = 1 & (self.super_states[index]>>shift);
                if bit > 0 { indexes.push((index*64)+shift) }
            }
        }

        if indexes.len() == 1 {
            self.collapsed = true;
            return;
        }
        //println!("collapse from {:?}",self.super_states);
        //println!("index from {:?}",indexes);

        if indexes.len() == 0 {panic!("trying to collapse 0 at {} {}",self.x,self.y)}

        //println!("pre collapse len {}",indexes.len());

        let mut r = rand::thread_rng();
        for _ in 0..(indexes.len()-1) {
            let i = r.gen_range(0..indexes.len());

            indexes.remove(i);
        }

        //println!("collapsed to {}",indexes[0]);

        let mut new_states = vec![0;self.super_states.len()];

        let byte_i = indexes[0]/64;
        let bit_i = indexes[0]%64;

        new_states[byte_i] = 1 << bit_i;

        self.super_states = new_states;

        //println!("post collapse {}",indexes.len());

        //println!("collapse at {} {} ,{}",self.x,self.y,self.super_states[0].index);
        self.collapsed = true;
        //println!("c_l_n:{}",self.super_states.len());

    }

    fn state(&self)->Rgba<u8>{
        for index in 0..self.super_states.len() {
            let byte = self.super_states[index];
            if byte == 0 {continue}
            for shift in 0..64 {
                let bit = 1 & (byte >> shift);
                if bit > 0 {
                    let index = (index * 64) + shift;
                    return self.lib[index].state[4];
                }
            }
        }
        return self.lib[0].state[4];
    }

    fn possible_states(&self, off : usize) -> Vec<usize>{

        let off = off;

        let mut lib_filter = Vec::new();
        //println!("lib {}",self.lib.len());
        //println!("s_s {:?}",self.super_states);
        for (s_i,u) in self.super_states.iter().enumerate() {
            //println!("bit {} {}",j,u);

            for b in 0..64 {
                let bit = u & (1<<b);
                if bit > 0 {
                    let l_i = (64*s_i)+b;
                    //println!("{}",l_i);

                    lib_filter.push(&self.lib[l_i]);
                }
            }
        }

        //println!("lib {}",lib_filter.len());

        //println!("states {}",self.super_states.len());
        
        let poss_states = lib_filter.iter().fold(vec![0;self.super_states.len()],|s_1,s_2|{
            let mut out = Vec::new();
            //println!("s_1 {}",s_1.len());
            for k in 0..s_1.len() {
                //println!("n_l {}",s_2.neighbours[off].len());
                let union = s_1[k] | s_2.neighbours[off][k];

                out.push( union);
            }

            out

        });

        //println!("poss {}",poss_states.len());

        poss_states
    }

    fn propograte(&mut self, super_states: Vec<usize>) -> bool{

        if self.super_states.len() == 0 {panic!("No States??");}

        if self.collapsed {return false;}

        let old_len = self.super_states.clone();

        for j in 0..super_states.len() {
            self.super_states[j] &= super_states[j];
        }

        let new_len = self.super_states.clone();

        if new_len == old_len {return false}

        true

    }
}


