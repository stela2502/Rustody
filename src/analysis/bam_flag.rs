

#[derive(Debug, Clone )]
pub struct BamFlag {
    pub paired: bool,              // 0x1
    pub proper_pair: bool,         // 0x2
    pub unmapped: bool,            // 0x4
    pub mate_unmapped: bool,       // 0x8
    pub reverse_strand: bool,      // 0x10
    pub mate_reverse_strand: bool, // 0x20
    pub read1: bool,               // 0x40
    pub read2: bool,               // 0x80
    pub secondary: bool,           // 0x100 - multimapper
    pub qc_fail: bool,             // 0x200
    pub duplicate: bool,           // 0x400
}

impl Default for BamFlag {
    fn default() -> Self {
        Self {
            paired: false,
            proper_pair: false,
            unmapped: false,
            mate_unmapped: false,
            reverse_strand: false,
            mate_reverse_strand: false,
            read1: false,
            read2: false,
            secondary: false,
            qc_fail: false,
            duplicate: false,
        }
    }
}

impl BamFlag {

    /// create a flag based on its numeric value
    pub fn new(flag: u16) -> Self {
        Self {
            paired: flag & 0x1 != 0,
            proper_pair: flag & 0x2 != 0,
            unmapped: flag & 0x4 != 0,
            mate_unmapped: flag & 0x8 != 0,
            reverse_strand: flag & 0x10 != 0,
            mate_reverse_strand: flag & 0x20 != 0,
            read1: flag & 0x40 != 0,
            read2: flag & 0x80 != 0,
            secondary: flag & 0x100 != 0,
            qc_fail: flag & 0x200 != 0,
            duplicate: flag & 0x400 != 0,
        }
    }    


    /// Convert the boolean flags into the corresponding SAM flag integer
    pub fn to_sam(&self) -> u16 {
        let mut flag = 0;
        if self.paired { flag |= 0x1; }
        if self.proper_pair { flag |= 0x2; }
        if self.unmapped { flag |= 0x4; }
        if self.mate_unmapped { flag |= 0x8; }
        if self.reverse_strand { flag |= 0x10; }
        if self.mate_reverse_strand { flag |= 0x20; }
        if self.read1 { flag |= 0x40; }
        if self.read2 { flag |= 0x80; }
        if self.secondary { flag |= 0x100; }
        if self.qc_fail { flag |= 0x200; }
        if self.duplicate { flag |= 0x400; }
        flag
    }

    /// simple way to check the deffernet values like flag.is("paired")
    pub fn is(&self, tag: &str) -> bool {
        match tag {
            "paired" => self.paired,
            "proper_pair" => self.proper_pair,
            "unmapped" => self.unmapped,
            "mate_unmapped" => self.mate_unmapped,
            "reverse_strand" => self.reverse_strand,
            "mate_reverse_strand" => self.mate_reverse_strand,
            "read1" => self.read1,
            "read2" => self.read2,
            "secondary" => self.secondary,
            "qc_fail" => self.qc_fail,
            "duplicate" => self.duplicate,
            _ => false, // Return false for unknown tags
        }
    }

    /// Set a specific flag
    pub fn set_paired(&mut self, value: bool) {
        self.paired = value;
    }
    pub fn set_proper_pair(&mut self, value: bool) {
        self.proper_pair = value;
    }
    pub fn set_unmapped(&mut self, value: bool) {
        self.unmapped = value;
    }
    pub fn set_mate_unmapped(&mut self, value: bool)  {
        self.mate_unmapped = value;
    }
    pub fn set_reverse_strand(&mut self, value: bool) {
        self.reverse_strand = value;
    }
    pub fn set_mate_reverse_strand(&mut self, value: bool) {
        self.mate_reverse_strand = value;
    }
    pub fn set_read1(&mut self, value: bool)  {
        self.read1 = value;
    }
    pub fn set_read2(&mut self, value: bool) {
        self.read2 = value;
    }
    pub fn set_secondary(&mut self, value: bool) {
        self.secondary = value;
    }
    pub fn set_qc_fail(&mut self, value: bool) {
        self.qc_fail = value;
    }
    pub fn set_duplicate(&mut self, value: bool)  {
        self.duplicate = value;
    }

}

