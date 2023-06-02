use std::cmp::min;
use std::iter::zip;
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

        library.push(State{state,neighbours:Default::default(),index:l_i});
        l_i += 1;
    }

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

        let min_x = min_index%width;
        let min_y = min_index/width;

        println!("e_i:{}",min_index);

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
                    if c.propograte(x_off, y_off, u_c.possible_states(i as usize).clone()) {
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
    state: [Rgba<u8>;9]
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
            let b_i = (other.index + 64 - 1)/64;
            let b_o = other.index%64;

            while self.neighbours[i].len() <= b_i {self.neighbours[i].push(0)}

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

        c.super_states[len-1] = (2 << (size%64))-1;

        c
    }

    fn collapse_coef(&self) -> usize{
        let s = self.super_states.iter().map(|u|u.count_ones()).reduce(|a,b|a+b);

        if self.collapsed {
            return usize::MAX;
        }

        s.unwrap() as usize
    }

    fn collapse(&mut self){
        //println!("c_l:{}",self.super_states.len());
        if self.super_states.len() == 1 {
            self.collapsed = true;
            return;
        }

        if self.super_states.len() == 0 {panic!("trying to collapse 0 at {} {}",self.x,self.y)}

        let mut r = rand::thread_rng();
        for _ in 0..(self.super_states.len()-1) {
            let i = r.gen_range(0..self.super_states.len());
            self.super_states.remove(i);
        }
        //println!("collapse at {} {} ,{}",self.x,self.y,self.super_states[0].index);
        self.collapsed = true;
        //println!("c_l_n:{}",self.super_states.len());

    }

    fn state(&self)->Rgba<u8>{
        let index = self.super_states[0];
        self.lib[index].state[4]
    }

    fn possible_states(&self, i : usize) -> Vec<usize>{

        let mut lib_filter = Vec::new();
        println!("lib {}",self.lib.len());
        println!("s_s {:?}",self.super_states);
        for (j,u) in self.super_states.iter().enumerate() {
            println!("bit {} {}",j,u);

            for b in 0..64 {
                let bit = u & (1<<b);
                if(bit > 0){
                    let l_i = (64*j)+b;
                    println!("{}",l_i);

                    lib_filter.push(self.lib[l_i].clone());
                }
            }
        }
        
        let poss_states = lib_filter.iter().fold(vec![0,lib_filter.len()],|s_1,s_2|{
            let mut out = Vec::new();
            for k in 0..s_1.len() {
                out.push( s_1[k] | s_2.neighbours[i][k]);
            }

            out

        });

        poss_states
    }

    fn propograte(&mut self,x_off:i32,y_off:i32,super_states: Vec<usize>) -> bool{

        if self.super_states.len() == 0 {panic!("No States??");}

        //let states : Vec<Rgba<u8>> = super_states.iter().map(|&s|s[4].clone()).collect();

        if self.collapsed {return false;}
        let i = x_off + (y_off * 3) + 4;

        let old_len = self.super_states.clone();

        //self.super_states = self.super_states.iter().filter(|s_s|super_states.contains(&s_s.index)).map(|s|s.clone()).collect();

        for j in 0..super_states.len() {
            self.super_states[j] &= super_states[j];
        }

        /*
        let new_len = self.super_states.len();
        if new_len == 0{
            println!("x:{} y:{}",self.x,self.y);
            for j in 0..196 as usize {
                if j % 14 == 0
                {
                    println!("");
                    print!("# ");
                }

                let states : Vec<usize> = old_len.iter().map(|s|s.index).collect();

                if super_states.contains(&j) {print!("o")}
                else if states.contains(&j) {print!("x")}
                else { print!(".") }
            }
            println!("");
            panic!("What? {:?}",super_states);
        }

        if new_len == old_len.len() {return false;}




        println!("x:{} y:{}",self.x,self.y);
        for j in 0..196 as usize {
            if j % 14 == 0
            {
                println!("");
                print!("# ");
            }

            let states : Vec<usize> = old_len.iter().map(|s|s.index).collect();

            if super_states.contains(&j) {print!("o")}
            else if states.contains(&j) {print!("x")}
            else { print!(".") }
        }
        println!("");

         */

        //let self_states : Vec<Rgba<u8>> = super_states.iter().map(|&s|s[i as usize].clone()).collect();

        //println!("i:{} x:{} y:{}",i,x_off,y_off);
        /*
        let old_states = self.super_states.clone();
        let squares : Vec<(usize,(usize,Rgba<u8>))> = self.super_states.iter().enumerate().map(|(j,&(k,s))| (j,(k,s[i as usize]))).collect();
        let squares_filter : Vec<&usize> = squares.iter().filter(|&&(_,(_,s))| states.contains(&s)).map(|(i,_)| i).collect();

        let self_square : Vec<(usize,(usize,Rgba<u8>))> = self.super_states.iter().enumerate().map(|(j,&(k,s))| (j,(k,s[4]))).collect();
        let mut self_filter : Vec<&usize> = squares.iter().filter(|&&(_,(_,s))| self_states.contains(&s)).map(|(i,_)| i).collect();

        self_filter.extend(squares_filter);

        let squares_filter : HashSet<&usize> = self_filter.iter().map(|&u|u).clone().collect();


        if squares_filter.len() == old_states.len() {return false;}
        //println!("x:{} y{}",self.x,self.y);
        //println!("old: l:{} {:?}",old_states.len(),old_states.iter().map(|(i,_)|i).collect::<Vec<&usize>>());
        //println!("other: l_s:{} l_f:{}",squares.len(),squares_filter.len());

        //println!("s:{} s_f:{} s_t:{}",squares.len(),squares_filter.len(),states.len());
        //println!("i: {} : {:?}",i,states);
        //println!("squares:{:?}\nstates:{:?}\noff:{},{}\nfilter:{:?}",self.super_states,states,x_off,y_off,squares_filter);

        let indexes : Vec<(usize,&(usize,[Rgba<u8>;9]))> = old_states.iter().enumerate().collect();
        let removed : Vec<(&usize,&usize)> = indexes.iter().filter(|(i,_)|!squares_filter.contains(&i)).map(|(i,(r,_))|(i,r)).collect();
        //let keep : Vec<(&usize,&usize)> = indexes.iter().filter(|(i,_)|squares_filter.contains(&i)).map(|(i,(r,_))|(i,r)).collect();


        for j in 0..196 {
            if j % 14 == 0
            {
                println!("");
                print!("# ");
            }
            let r_i : Vec<&(&usize,&usize)> = removed.iter().filter(|(_,&k)|j == k).collect();
            let s_i : Vec<&(&usize,&usize)> = keep.iter().filter(|(_,&k)| j == k).collect();
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

        //println!("n:{}",self.super_states.len());
        //println!("c:{}",changed);
*/
        true

    }
}


