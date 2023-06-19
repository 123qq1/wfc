use std::cmp::min;
use std::collections::HashMap;
use image::{DynamicImage, GenericImageView, Rgba, GenericImage};
use image::io::Reader as ImageReader;
use rand::distributions::WeightedIndex;
use rand::prelude::*;

fn main() {
    println!("Hello, wave function collapse!");

    let img = ImageReader::open("img/OfficeBuilding.png").unwrap().decode().unwrap();

    let mut library = create_lib(img);
    let mut poss_state_cache : HashMap<Vec<usize>,[Vec<usize>;9]> = HashMap::new();

    println!("Sampling done {}",library.len());

    let lib_c = library.clone();

    //fill out the library with neighbours
    for state in &mut library {
        for neighbour in &lib_c {
            for i in 0..9 {
                if i == 4 {continue;}
                state.neighbour_check(neighbour,i);
            }
        }
    }

    for state in &mut library {
        state.find_stability();
    }

    let library = library;

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
    let height = 50;
    let width = 50;
    let size = (height * width) as usize;

    let mut out = vec![Vec::with_capacity(width);height];

    for x in 0..width as u32 {
        for y in 0..height as u32 {
            out[y as usize].push(Cell::new(&library,x,y,library.len()));
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

    let mut first = true;

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
        let (mut min_index,_) = entropys.iter().enumerate().min_by_key(|(_,&t)|t).unwrap();

        if first {
            first = false;
            min_index = (width/2) + ((height/2)*64)
        }

        let min_x = min_index%width;
        let min_y = min_index/width;

        //println!("collapsing {} {}",min_x,min_y);

        out[min_y][min_x].collapse();
        out[min_y][min_x].update_poss_states(&mut poss_state_cache);
        updates.push((min_y,min_x));

        while updates.len() > 0 {
            let _pre_len = updates.len();

            let (u_y, u_x) = updates.pop().unwrap();

            //println!("working on: {} {}",u_x,u_y);

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

                    //let x_c_off = c.x as i32 - min_x as i32;
                    //let y_c_off = c.y as i32 - min_y as i32;

                    //if x_c_off.abs() > 1 {continue;}
                    //if y_c_off.abs() > 1 {continue;}

                    //println!("x:{:02} y:{:02}",x_off,y_off);

                    let _n_s = (y_2 - y_1) * (x_2 - x_1);
                    //println!("{}",n_s);
                    //println!("{} {}",c.y,c.x);
                    if c.propograte( u_c.possible_states(i as usize)) {
                        c.update_poss_states(&mut poss_state_cache);
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
    println!("done");

    let mut d_image = DynamicImage::new_rgba8(width as u32,height as u32);

    for w in 0..width {
        for h in 0..height {

            d_image.put_pixel(w as u32,h as u32,out[h][w].state());

        }
    }
    let p = "./img/out.png";
    d_image.save(p).unwrap();
}

fn create_lib(img: DynamicImage) -> Vec<State> {
    let state_size = ((img.height() - 2) * (img.width() - 2)) as usize;

    let mut library = Vec::with_capacity(state_size);

    for pixel in img.pixels() {
        if pixel.0 == 0 || pixel.1 == 0 { continue; }
        if pixel.0 == (img.width() - 1) || pixel.1 == (img.height() - 1) { continue; }

        let x = pixel.0;
        let y = pixel.1;

        let state = [
            img.get_pixel(x - 1, y - 1), img.get_pixel(x, y - 1), img.get_pixel(x + 1, y - 1),
            img.get_pixel(x - 1, y), img.get_pixel(x, y), img.get_pixel(x + 1, y),
            img.get_pixel(x - 1, y + 1), img.get_pixel(x, y + 1), img.get_pixel(x + 1, y + 1)
        ];
        if !library.contains(&state) {
            library.push(state);
        }
    }

    let state_size = library.len();
    let req_bytes = (state_size + 64 - 1) / 64;

    let n_base =
        [
            vec![0; req_bytes], vec![0; req_bytes], vec![0; req_bytes],
            vec![0; req_bytes], vec![0; req_bytes], vec![0; req_bytes],
            vec![0; req_bytes], vec![0; req_bytes], vec![0; req_bytes]
        ];

    let library: Vec<State> = library.iter().enumerate().map(|(i, s)|
        State { neighbours: n_base.clone(), state: s.clone(), index: i, /*size: state_size*/stability:0 }
    ).collect();
    library
}

#[derive(Clone,Debug)]
struct State{
    index : usize,
    neighbours: [Vec<usize>;9],
    state: [Rgba<u8>;9],
    stability: usize,
    //size: usize
}

impl PartialEq for State{
    fn eq(&self, other: &Self) -> bool {
        self.state == other.state
    }
}

impl State{

    fn find_stability(&mut self){
        let not_self = self.neighbours.iter().enumerate().filter(|(i,_)|*i != 4).map(|(_,n)|n);

        self.stability = not_self.fold(usize::MAX,|n_1,n_2|{
            let ones = n_2.iter().fold(0,|b_1,b_2|{b_1 + b_2.count_ones()}) as usize;
            min(n_1,ones)
            //n_1+ones
        });
        //println!("s {}",self.stability);
    }

    fn neighbour_check(&mut self, other :&State, i : usize){
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
struct Cell<'a>{
    lib: &'a Vec<State>,
    super_states: Vec<usize>,
    poss_states: [Vec<usize>;9],
    x: u32,
    y: u32,
    collapsed: bool,
    //size: usize,
}

impl Cell<'_>{
    fn new(lib : &Vec<State>,x : u32,y: u32,size: usize) -> Cell{

        let len = (size + 64 - 1)/64;

        let poss_states =
            [
                vec![0;len],vec![0;len],vec![0;len],
                vec![0;len],vec![0;len],vec![0;len],
                vec![0;len],vec![0;len],vec![0;len],
            ];

        let mut c = Cell{lib,x,y,collapsed:false,/*size,*/super_states:vec![usize::MAX;len],poss_states};

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
        let mut weights = Vec::new();

        for index in 0..self.super_states.len() {
            for shift in 0..64 {
                let bit = 1 & (self.super_states[index]>>shift);
                if bit > 0 {
                    let s_index = (index*64)+shift;
                    indexes.push(s_index);
                    weights.push(self.lib[s_index].stability)
                }
            }
        }

        if indexes.len() <= 1 {
            //self.update_poss_states();
            self.collapsed = true;
            return;
        }


        //println!("pre collapse len {}",indexes.len());
        let dist = WeightedIndex::new(&weights).unwrap();
        let mut r = rand::thread_rng();

        let chosen = dist.sample(&mut r);

        //println!("collapse stability {}",weights[chosen]);

        //println!("collapsed to {}",indexes[0]);

        let mut new_states = vec![0;self.super_states.len()];

        let byte_i = indexes[chosen]/64;
        let bit_i = indexes[chosen]%64;

        new_states[byte_i] = 1 << bit_i;

        self.super_states = new_states;

        //println!("post collapse {}",indexes.len());

        //println!("collapse at {} {} ,{}",self.x,self.y,self.super_states[0].index);
        self.collapsed = true;
        //self.update_poss_states();
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

    fn update_poss_states(&mut self, cache : &mut HashMap<Vec<usize>,[Vec<usize>;9]>){

        if let Some(poss_state) = cache.get::<Vec<usize>>(&self.super_states){
            self.poss_states = poss_state.clone();
            return;
        }


        for n_c in 0..9 {
            for b_i in 0..self.super_states.len() {
                self.poss_states[n_c][b_i] = 0;
            }
        }

        for (s_i,byte) in self.super_states.iter().enumerate() {
            for b in 0..64 {
                let bit = byte & (1<<b);
                if bit > 0 {
                    let l_i = (64*s_i)+b;
                    for n_i in 0..9 {
                        for (n_b_i, n_byte) in self.lib[l_i].neighbours[n_i].iter().enumerate() {
                            self.poss_states[n_i][n_b_i] |= n_byte;
                        }
                    }
                }
            }
        }

        if !cache.contains_key::<Vec<usize>>(&self.super_states) {
            cache.insert(self.super_states.clone(),self.poss_states.clone());
        }
    }

    fn possible_states(&self, off : usize) -> &Vec<usize>{


        return &self.poss_states[off];
        /*
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

                    lib_filter.push(&self.lib[l_i].neighbours[off]);
                }
            }
        }

        //println!("lib {}",lib_filter.len());

        //println!("states {}",self.super_states.len());
        
        let mut poss_states = vec![0;self.super_states.len()];

        for s in lib_filter {
            for (i,n) in s.iter().enumerate() {
                poss_states[i] |= n;
            }
        }

        //println!("poss {}",poss_states.len());

        &poss_states

         */
    }

    fn propograte(&mut self, super_states: &Vec<usize>) -> bool{

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


