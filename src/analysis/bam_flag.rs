

#[derive(Debug, Default)]
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

impl BamFlag {
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

    /// Set a flag based on its numeric value
    pub fn set(&mut self, flag: u16) {
        if flag & 0x1 != 0 { self.paired = true; }
        if flag & 0x2 != 0 { self.proper_pair = true; }
        if flag & 0x4 != 0 { self.unmapped = true; }
        if flag & 0x8 != 0 { self.mate_unmapped = true; }
        if flag & 0x10 != 0 { self.reverse_strand = true; }
        if flag & 0x20 != 0 { self.mate_reverse_strand = true; }
        if flag & 0x40 != 0 { self.read1 = true; }
        if flag & 0x80 != 0 { self.read2 = true; }
        if flag & 0x100 != 0 { self.secondary = true; }
        if flag & 0x200 != 0 { self.qc_fail = true; }
        if flag & 0x400 != 0 { self.duplicate = true; }
    }

}

